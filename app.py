
from flask import Flask, render_template, request, send_file

import pandas as pd
import numpy as np
import requests
import os

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

import matplotlib.pyplot as plt
import gseapy as gp

from adjustText import adjust_text
from bs4 import BeautifulSoup

app = Flask(__name__)


# =========================================================
# FILE LOADER
# =========================================================

def load_file(file):

    filename = file.filename.lower()

    # ================= EXCEL =================

    if filename.endswith(".xlsx") or filename.endswith(".xls"):

        try:

            df = pd.read_excel(file)

            df = df.dropna(axis=1, how="all")

            return df

        except Exception as e:

            raise ValueError(
                f"Excel read error: {e}"
            )

    # ================= CSV / TSV / TXT =================

    separators = [",", "\t", ";"]

    for sep in separators:

        try:

            file.seek(0)

            df = pd.read_csv(

                file,

                sep=sep,

                engine="python",

                comment="#",

                on_bad_lines="skip",

                compression="infer"

            )

            df = df.dropna(axis=1, how="all")

            if df.shape[1] >= 2:

                return df

        except:
            continue

    raise ValueError(

        "Could not read file. "

        "Please upload CSV, TSV, TXT, XLSX, or GZ files."

    )


# =========================================================
# FETCH GENE INFO
# =========================================================

def fetch_gene_info(gene):

    try:

        url = "https://mygene.info/v3/query"

        params = {

            "q": f"symbol:{gene}",

            "species": "human",

            "fields": "symbol,name,summary",

            "size": 1

        }

        response = requests.get(url, params=params)

        data = response.json()

        if data.get("hits"):

            hit = data["hits"][0]

            return (

                hit.get("symbol", gene),

                hit.get("name", "Unknown"),

                hit.get("summary", "No summary available")

            )

    except:
        pass

    return gene, "Unknown", "No summary available"


# =========================================================
# DRUGGABILITY
# =========================================================

def fetch_drug_info(gene, summary):

    try:

        url = f"https://dgidb.org/api/v2/interactions.json?genes={gene}"

        r = requests.get(url, timeout=3)

        if r.status_code == 200:

            data = r.json()

            matches = data.get("matchedTerms", [])

            if matches:

                interactions = matches[0].get("interactions", [])

                if interactions:

                    return "Yes"

    except:
        pass

    text = str(summary).lower()

    if any(k in text for k in [

        "receptor",

        "kinase",

        "enzyme",

        "channel",

        "gpcr"

    ]):

        return "Yes"

    if any(k in text for k in [

        "tubulin",

        "microtubule",

        "cytoskeleton"

    ]):

        return "Potential"

    return "No"


# =========================================================
# PRIORITIZATION
# =========================================================

def generate_reason(log2fc, druggable, summary):

    reasons = []

    if log2fc > 2:

        reasons.append("Strong overexpression")

    elif log2fc > 1:

        reasons.append("Moderate overexpression")

    if druggable == "Yes":

        reasons.append("High druggability potential")

    elif druggable == "Potential":

        reasons.append("Potential therapeutic relevance")

    text = str(summary).lower()

    if "receptor" in text:

        reasons.append("Receptor-associated target")

    if "kinase" in text:

        reasons.append("Kinase signaling involvement")

    if "enzyme" in text:

        reasons.append("Enzymatic target class")

    if "microtubule" in text or "tubulin" in text:

        reasons.append("Cytoskeletal involvement")

    if "immune" in text or "interleukin" in text:

        reasons.append("Immune-related pathway")

    if len(reasons) == 0:

        reasons.append("Transcriptomic significance")

    return " + ".join(reasons)


# =========================================================
# PUBMED
# =========================================================

def generate_pubmed_link(gene):

    query = str(gene).replace(" ", "+")

    return (

        f"https://pubmed.ncbi.nlm.nih.gov/?term={query}"

    )


def get_pubmed_evidence(gene):

    try:

        url = (

            f"https://pubmed.ncbi.nlm.nih.gov/"

            f"?term={gene}"

        )

        response = requests.get(

            url,

            timeout=5

        )

        soup = BeautifulSoup(

            response.text,

            "html.parser"

        )

        count_tag = soup.find(

            "meta",

            attrs={"name": "log_resultcount"}

        )

        if count_tag:

            count = int(

                count_tag["content"]

            )

        else:

            count = 0

    except:

        count = 0

    if count > 5000:

        return "Highly studied gene"

    elif count > 500:

        return "Moderately studied gene"

    else:

        return "Understudied target"


# =========================================================
# PATHWAY ANALYSIS
# =========================================================

def pathway_enrichment(df):

    genes = df[

        df["adj_p_value"] < 0.1

    ].index.tolist()

    if len(genes) < 5:

        return pd.DataFrame()

    try:

        enr = gp.enrichr(

            gene_list=genes,

            gene_sets=["KEGG_2021_Human"],

            organism="human",

            outdir=None

        )

        return enr.results.head(5)

    except:

        return pd.DataFrame()


# =========================================================
# MAIN ANALYSIS
# =========================================================

def analyze(file):

    df = load_file(file)

    # =====================================================
    # FIND GENE COLUMN
    # =====================================================

    gene_col = None

    for col in df.columns:

        sample_values = (

            df[col]

            .dropna()

            .astype(str)

            .head(20)

        )

        non_numeric = 0

        for val in sample_values:

            try:
                float(val)

            except:
                non_numeric += 1

        if len(sample_values) > 0:

            ratio = (

                non_numeric /

                len(sample_values)

            )

            if ratio > 0.8:

                gene_col = col

                break

    if gene_col is None:

        raise ValueError(

            "No gene column found."

        )

    # =====================================================
    # SPECIES VALIDATION
    # =====================================================

    gene_values = (

        df[gene_col]

        .head(20)

        .tolist()

    )

    gene_values = [

        str(g)

        for g in gene_values

        if pd.notna(g)

    ]

    # ================= FLY =================

    if any(

        g.startswith("FBgn")

        for g in gene_values

    ):

        raise ValueError(

            "Fly dataset detected. "

            "Current version supports human data only."

        )

    # ================= MOUSE =================

    if any(

        g.startswith("ENSMUSG")

        for g in gene_values

    ):

        raise ValueError(

            "Mouse dataset detected. "

            "Current version supports human data only."

        )

    # =====================================================
    # SET INDEX
    # =====================================================

    df.set_index(gene_col, inplace=True)

    # Remove duplicate genes

    df = df[~df.index.duplicated(keep="first")]

    # =====================================================
    # KEEP NUMERIC
    # =====================================================

    numeric_df = df.select_dtypes(

        include=[np.number]

    )

    if numeric_df.shape[1] < 2:

        raise ValueError(

            "Not enough numeric expression columns found."

        )

    df = numeric_df

    # =====================================================
    # REMOVE EMPTY
    # =====================================================

    df.dropna(axis=1, how="all", inplace=True)

    # =====================================================
    # SPLIT GROUPS
    # =====================================================

    cols = list(df.columns)

    mid = len(cols) // 2

    A = cols[:mid]

    B = cols[mid:]

    # =====================================================
    # MEANS
    # =====================================================

    df["Mean_A"] = df[A].mean(axis=1)

    df["Mean_B"] = df[B].mean(axis=1)

    # =====================================================
    # LOG2FC
    # =====================================================

    df["log2FC"] = np.log2(

        (df["Mean_B"] + 1) /

        (df["Mean_A"] + 1)

    )

    # =====================================================
    # P VALUES
    # =====================================================

    pvals = []

    for i in range(len(df)):

        try:

            a = df.iloc[i][A]

            b = df.iloc[i][B]

            _, p = ttest_ind(

                a,

                b,

                equal_var=False,

                nan_policy="omit"

            )

            if np.isnan(p):

                p = 1.0

        except:

            p = 1.0

        pvals.append(p)

    df["p_value"] = pvals

    # =====================================================
    # ADJUSTED P VALUE
    # =====================================================

    df["adj_p_value"] = multipletests(

        df["p_value"],

        method="fdr_bh"

    )[1]

    # =====================================================
    # VOLCANO
    # =====================================================

    df["neg_log10_p"] = -np.log10(

        df["p_value"] + 1e-10

    )

    df["category"] = "Not Significant"

    df.loc[

        (df["log2FC"] > 1) &

        (df["p_value"] < 0.05),

        "category"

    ] = "Upregulated"

    df.loc[

        (df["log2FC"] < -1) &

        (df["p_value"] < 0.05),

        "category"

    ] = "Downregulated"

    # =====================================================
    # FILTER
    # =====================================================

    filtered = df[

        (df["log2FC"].abs() > 1) &

        (df["p_value"] < 0.05)

    ]

    top = filtered.sort_values(

        by="log2FC",

        ascending=False

    ).head(10)

    if top.empty:

        raise ValueError(

            "No significant targets detected."

        )

    # =====================================================
    # ANNOTATION
    # =====================================================

    names = []
    infos = []
    reasons = []
    literature_links = []

    for gene in top.index:

        pname, pinfo, summary = fetch_gene_info(gene)

        dlabel = fetch_drug_info(

            gene,

            summary

        )

        reason = generate_reason(

            float(top.loc[gene, "log2FC"]),

            dlabel,

            summary

        )

        evidence = get_pubmed_evidence(gene)

        reason = (

            reason +

            " + " +

            evidence

        )

        pubmed = generate_pubmed_link(gene)

        names.append(pname)

        infos.append(pinfo)

        reasons.append(reason)

        literature_links.append(pubmed)

    top["Protein_Name"] = names
    top["Protein_Info"] = infos
    top["Why_Prioritized"] = reasons
    top["Literature"] = literature_links

    # =====================================================
    # TARGET SCORE
    # =====================================================

    top["log2FC_z"] = (

        top["log2FC"] -

        top["log2FC"].mean()

    ) / (

        top["log2FC"].std() + 1e-10

    )

    top["MeanB_z"] = (

        top["Mean_B"] -

        top["Mean_B"].mean()

    ) / (

        top["Mean_B"].std() + 1e-10

    )

    raw_score = (

        0.7 * top["log2FC_z"] +

        0.3 * top["MeanB_z"]

    )

    top["Target_Score"] = (

        (raw_score - raw_score.min()) /

        (raw_score.max() - raw_score.min() + 1e-10)

    )

    top = top.sort_values(

        by="Target_Score",

        ascending=False

    )

    top["Rank"] = range(

        1,

        len(top) + 1

    )

    # =====================================================
    # STATIC FOLDER
    # =====================================================

    if not os.path.exists("static"):

        os.makedirs("static")

    # =====================================================
    # SAVE EXCEL
    # =====================================================

    top.to_excel(

        "static/results.xlsx",

        index=False

    )

    # =====================================================
    # VOLCANO PLOT
    # =====================================================

    plt.figure(figsize=(8, 5))

    colors = {

        "Upregulated": "red",

        "Downregulated": "green",

        "Not Significant": "gray"

    }

    for cat in colors:

        subset = df[

            df["category"] == cat

        ]

        plt.scatter(

            subset["log2FC"],

            subset["neg_log10_p"],

            c=colors[cat],

            label=cat,

            alpha=0.7

        )

    plt.axvline(x=1, linestyle="--")
    plt.axvline(x=-1, linestyle="--")
    plt.axhline(y=-np.log10(0.05), linestyle="--")

    texts = []

    for gene in top.head(5).index:

        x = df.loc[gene, "log2FC"]

        y = df.loc[gene, "neg_log10_p"]

        texts.append(

            plt.text(

                x,

                y,

                gene,

                fontsize=8

            )

        )

    adjust_text(texts)

    plt.xlabel("log2 Fold Change")

    plt.ylabel("-log10(p-value)")

    plt.title("Volcano Plot")

    plt.legend()

    plt.tight_layout()

    volcano_path = "static/volcano.png"

    plt.savefig(volcano_path)

    plt.close()

    # =====================================================
    # PATHWAY PLOT
    # =====================================================

    pathway_df = pathway_enrichment(df)

    pathway_plot = None

    if not pathway_df.empty:

        plt.figure(figsize=(7, 4))

        plt.barh(

            pathway_df["Term"],

            -np.log10(

                pathway_df[

                    "Adjusted P-value"

                ]

            )

        )

        plt.xlabel("-log10(adj p-value)")

        plt.title("Top Enriched Pathways")

        plt.tight_layout()

        pathway_plot = "static/pathway.png"

        plt.savefig(pathway_plot)

        plt.close()

    # =====================================================
    # FINAL TABLE
    # =====================================================

    top = top[[

        "Rank",

        "Protein_Name",

        "Protein_Info",

        "Why_Prioritized",

        "Target_Score",

        "log2FC",

        "adj_p_value",

        "Literature"

    ]].round(3)

    return (

        top,

        volcano_path,

        pathway_plot

    )


# =========================================================
# MAIN ROUTE
# =========================================================

@app.route("/", methods=["GET", "POST"])

def index():

    if request.method == "POST":

        file = request.files["file"]

        result, plot, pathway_plot = analyze(file)

        result["Literature"] = result["Literature"].apply(

            lambda x:
            f'<a href="{x}" target="_blank">PubMed</a>'

        )

        styled = result.style.background_gradient(

            subset=["Target_Score"],

            cmap="viridis"

        )

        return render_template(

            "index.html",

            table=styled.to_html(escape=False),

            plot=plot,

            pathway_plot=pathway_plot

        )

    return render_template("index.html")


# =========================================================
# DOWNLOAD RESULTS
# =========================================================

@app.route("/download")

def download_results():

    return send_file(

        "static/results.xlsx",

        as_attachment=True

    )


# =========================================================
# RUN
# =========================================================

if __name__ == "__main__":

    app.run(debug=True)
