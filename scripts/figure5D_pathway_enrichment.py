import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
os.makedirs(FIG5, exist_ok=True)

de_file = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
de = pd.read_csv(de_file)

# ------------------------------------------------------------
# Select genes
# ------------------------------------------------------------
meta_genes = (
    de[(de["direction"] == "metastasis_up")]
    .sort_values("log2FC", ascending=False)["gene"]
    .head(300)
    .tolist()
)

primary_genes = (
    de[(de["direction"] == "primary_up")]
    .sort_values("log2FC", ascending=True)["gene"]
    .head(300)
    .tolist()
)

# CLEAN genes
meta_genes = [g for g in meta_genes if isinstance(g, str)]
primary_genes = [g for g in primary_genes if isinstance(g, str)]

print("Metastasis genes:", len(meta_genes))
print("Primary genes:", len(primary_genes))

# ------------------------------------------------------------
# Gene sets
# ------------------------------------------------------------
gene_sets = [
    "GO_Biological_Process_2023",
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human"
]

# ------------------------------------------------------------
# Enrichment function
# ------------------------------------------------------------
def run_enrichment(genes, label):
    all_results = []

    for gs in gene_sets:
        print(f"Running {label} | {gs}")

        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gs,
                organism="human",  # FIXED
                outdir=None,
                cutoff=1.0
            )

            res = enr.results.copy()

            if res is None or res.shape[0] == 0:
                print("No results for", gs)
                continue

            res["source_gene_set"] = gs
            res["comparison"] = label
            all_results.append(res)

        except Exception as e:
            print("ERROR in", gs, ":", e)

    if len(all_results) == 0:
        return None, None

    out = pd.concat(all_results, ignore_index=True)

    # Clean columns
    out["Adjusted P-value"] = pd.to_numeric(out["Adjusted P-value"], errors="coerce")
    out["Combined Score"] = pd.to_numeric(out["Combined Score"], errors="coerce")

    out = out.dropna(subset=["Adjusted P-value"])
    out = out.sort_values("Adjusted P-value")

    out_file = os.path.join(FIG5, f"Figure5D_{label}_all_pathways.csv")
    out.to_csv(out_file, index=False)

    top = out.head(20).copy()

    top_file = os.path.join(FIG5, f"Figure5D_{label}_top20.csv")
    top.to_csv(top_file, index=False)

    return out, top

# ------------------------------------------------------------
# Run enrichment
# ------------------------------------------------------------
meta_all, meta_top = run_enrichment(meta_genes, "metastasis_up")
primary_all, primary_top = run_enrichment(primary_genes, "primary_up")

# ------------------------------------------------------------
# Plot function
# ------------------------------------------------------------
def plot_bar(top_df, label, color, filename):
    if top_df is None or top_df.shape[0] == 0:
        print("No data to plot for", label)
        return

    df = top_df.copy()

    df["minus_log10_adj_p"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-300))
    df = df.sort_values("minus_log10_adj_p", ascending=True)

    plt.figure(figsize=(10, 7))

    plt.barh(df["Term"], df["minus_log10_adj_p"], color=color)

    plt.xlabel("-log10 adjusted p-value")
    plt.ylabel("")
    plt.title(label.replace("_", " ") + " pathway enrichment")

    plt.tight_layout()

    plt.savefig(os.path.join(FIG5, filename + ".tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG5, filename + ".pdf"), bbox_inches="tight")

    plt.close()

# ------------------------------------------------------------
# Plot both
# ------------------------------------------------------------
plot_bar(
    meta_top,
    "metastasis_up",
    "#D73027",
    "Figure5D_metastasis_barplot"
)

plot_bar(
    primary_top,
    "primary_up",
    "#4575B4",
    "Figure5D_primary_barplot"
)

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
report_path = os.path.join(FIG5, "Figure5D_PATHWAY_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 5D PATHWAY ENRICHMENT REPORT\n\n")

    f.write("Metastasis genes used:\n")
    f.write(", ".join(meta_genes[:100]) + "\n\n")

    f.write("Primary genes used:\n")
    f.write(", ".join(primary_genes[:100]) + "\n\n")

    if meta_top is not None:
        f.write("Top metastasis pathways:\n")
        f.write(meta_top[["Term", "Adjusted P-value", "Combined Score"]].to_string(index=False))
        f.write("\n\n")

    if primary_top is not None:
        f.write("Top primary pathways:\n")
        f.write(primary_top[["Term", "Adjusted P-value", "Combined Score"]].to_string(index=False))
        f.write("\n\n")

    f.write("Files:\n")
    for x in sorted(os.listdir(FIG5)):
        if "Figure5D" in x:
            f.write(x + "\n")

print("\nDONE FIGURE 5D")
print(report_path)

if os.path.exists(report_path):
    print("\n===== REPORT PREVIEW =====\n")
    print(open(report_path, encoding="utf-8").read())