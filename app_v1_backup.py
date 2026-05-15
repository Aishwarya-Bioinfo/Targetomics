from flask import Flask, render_template, request
import pandas as pd
import numpy as np

app = Flask(__name__)

# Protein mapping (basic version)
protein_info = {
    "FGF19": "Fibroblast growth factor 19 (Growth factor)",
    "TRPA1": "Ion channel receptor",
    "GDNF": "Neurotrophic factor",
    "DCC": "Tumor suppressor receptor",
    "MC4R": "G-protein coupled receptor",
    "CHGB": "Secretory protein",
    "ADAM2": "Metalloprotease",
    "LHX1": "Transcription regulator",
    "CARTPT": "Neuropeptide",
    "NEUROD1": "Transcription factor"
}

def analyze(file):
    # Load file
    df = pd.read_csv(file, sep="\t")

    # Set gene name as index
    df.set_index("GeneName", inplace=True)

    # Separate groups
    A_cols = [col for col in df.columns if col.startswith("A_")]
    B_cols = [col for col in df.columns if col.startswith("B_")]

    # Calculate means
    df["Mean_A"] = df[A_cols].mean(axis=1)
    df["Mean_B"] = df[B_cols].mean(axis=1)

    # Calculate log2FC
    df["log2FC"] = np.log2((df["Mean_B"] + 1) / (df["Mean_A"] + 1))

    # Filter
    filtered = df[(df["Mean_B"] > 50) & (df["log2FC"] > 2)]

    # Top genes
    top = filtered.sort_values(by="log2FC", ascending=False).head(10)

    # Add protein info
    top["Protein_Info"] = top.index.map(protein_info).fillna("Unknown protein")

    # Select columns
    top = top[["Protein_Info", "Mean_A", "Mean_B", "log2FC"]]

    return top


@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        file = request.files["file"]

        if file.filename == "":
            return "No file selected"

        result = analyze(file)

        return render_template("index.html", table=result.to_html())

    return render_template("index.html")


if __name__ == "__main__":
    app.run(debug=True)
