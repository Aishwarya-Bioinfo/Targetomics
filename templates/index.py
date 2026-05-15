<!DOCTYPE html>
<html>
<head>
    <title>Targetomics</title>

    <style>
        body {
            font-family: Arial, sans-serif;
            background-color: #f5f7fa;
            margin: 0;
        }

        .container {
            width: 85%;
            margin: 40px auto;
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0px 4px 12px rgba(0,0,0,0.1);
        }

        h1 {
            margin-bottom: 5px;
            color: #2c3e50;
        }

        h3 {
            color: #7f8c8d;
            font-weight: normal;
        }

        form {
            margin-top: 20px;
            margin-bottom: 20px;
        }

        button {
            background-color: #3498db;
            color: white;
            border: none;
            padding: 8px 18px;
            border-radius: 6px;
            cursor: pointer;
        }

        button:hover {
            background-color: #2980b9;
        }

        .download-btn {
            display: inline-block;
            margin-top: 15px;
            background-color: #27ae60;
            color: white;
            padding: 8px 16px;
            border-radius: 6px;
            text-decoration: none;
        }

        .download-btn:hover {
            background-color: #1e8449;
        }

        table {
            border-collapse: collapse;
            width: 100%;
            margin-top: 15px;
            font-size: 14px;
        }

        th {
            background-color: #2c3e50;
            color: white;
            padding: 10px;
        }

        td {
            padding: 8px;
            text-align: center;
        }

        tr:nth-child(even) {
            background-color: #f2f6fa;
        }

        tr:hover {
            background-color: #eaf2f8;
        }

        td:nth-child(2) {
            text-align: left;
            padding-left: 15px;
        }

        td:nth-child(3),
        td:nth-child(4),
        td:nth-child(5) {
            text-align: right;
            padding-right: 15px;
        }
    </style>
</head>

<body>

<div class="container">

    <h1>Targetomics</h1>
    <h3>Transcriptomics → Protein Target Prediction Tool</h3>

    <form method="POST" enctype="multipart/form-data">
        <input type="file" name="file" required>
        <button type="submit">Analyze</button>
    </form>

    {% if table %}
        <h2>Top Protein Targets</h2>
        {{ table|safe }}

        <a href="/download" class="download-btn">Download Results (CSV)</a>
    {% endif %}

</div>

</body>
</html>
