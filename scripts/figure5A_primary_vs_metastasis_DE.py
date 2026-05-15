import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
os.makedirs(FIG5, exist_ok=True)

adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")
adata = sc.read_h5ad(adata_path)

print("Loaded:")
print(adata)

# Keep only primary and metastasis cells
adata_pm = adata[adata.obs["label"].isin(["primary", "metastasis"])].copy()

print("Primary/metastasis subset:")
print(adata_pm)
print(pd.crosstab(adata_pm.obs["label"], adata_pm.obs["cell_type"]))

# Log-normalize for DE
adata_expr = adata_pm.copy()
sc.pp.normalize_total(adata_expr, target_sum=1e4)
sc.pp.log1p(adata_expr)

# Differential expression: metastasis vs primary
print("Running DE: metastasis vs primary...")

sc.tl.rank_genes_groups(
    adata_expr,
    groupby="label",
    groups=["metastasis"],
    reference="primary",
    method="wilcoxon",
    n_genes=adata_expr.n_vars
)

de = sc.get.rank_genes_groups_df(adata_expr, group="metastasis")

# Clean table
de = de.rename(columns={
    "names": "gene",
    "scores": "score",
    "logfoldchanges": "log2FC",
    "pvals": "pval",
    "pvals_adj": "padj"
})

de["neg_log10_padj"] = -np.log10(de["padj"].replace(0, np.nextafter(0, 1)))
de["direction"] = "not_significant"
de.loc[(de["padj"] < 0.05) & (de["log2FC"] > 0.25), "direction"] = "metastasis_up"
de.loc[(de["padj"] < 0.05) & (de["log2FC"] < -0.25), "direction"] = "primary_up"

de = de.sort_values(["padj", "log2FC"], ascending=[True, False])

de_out = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
de.to_csv(de_out, index=False)

# Top genes
top_meta = de[de["direction"] == "metastasis_up"].sort_values("log2FC", ascending=False).head(50)
top_primary = de[de["direction"] == "primary_up"].sort_values("log2FC", ascending=True).head(50)

top_meta.to_csv(os.path.join(FIG5, "Figure5A_top50_metastasis_up_genes.csv"), index=False)
top_primary.to_csv(os.path.join(FIG5, "Figure5A_top50_primary_up_genes.csv"), index=False)

# Volcano plot
plot_df = de.copy()
plot_df["plot_color"] = "gray"
plot_df.loc[plot_df["direction"] == "metastasis_up", "plot_color"] = "red"
plot_df.loc[plot_df["direction"] == "primary_up", "plot_color"] = "blue"

fig, ax = plt.subplots(figsize=(8, 7))

ax.scatter(
    plot_df["log2FC"],
    plot_df["neg_log10_padj"],
    c=plot_df["plot_color"],
    s=6,
    alpha=0.55,
    linewidths=0
)

ax.axvline(0.25, linestyle="--", linewidth=0.8)
ax.axvline(-0.25, linestyle="--", linewidth=0.8)
ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8)

ax.set_xlabel("log2 fold change: metastasis vs primary")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Primary vs metastasis differential expression")

# Label top genes
label_genes = pd.concat([
    top_meta.head(10),
    top_primary.head(10)
])

for _, row in label_genes.iterrows():
    ax.text(
        row["log2FC"],
        row["neg_log10_padj"],
        row["gene"],
        fontsize=7,
        ha="center",
        va="bottom"
    )

plt.tight_layout()

plt.savefig(
    os.path.join(FIG5, "Figure5A_primary_vs_metastasis_volcano.tiff"),
    dpi=300,
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(FIG5, "Figure5A_primary_vs_metastasis_volcano.pdf"),
    bbox_inches="tight"
)
plt.close()

# Summary report
report_path = os.path.join(FIG5, "Figure5A_DE_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 5A PRIMARY VS METASTASIS DIFFERENTIAL EXPRESSION REPORT\n\n")

    f.write("1. Subset size\n")
    f.write(str(adata_pm))
    f.write("\n\n")

    f.write("2. Label by cell type table\n")
    f.write(pd.crosstab(adata_pm.obs["label"], adata_pm.obs["cell_type"]).to_string())
    f.write("\n\n")

    f.write("3. DE gene counts\n")
    f.write(de["direction"].value_counts().to_string())
    f.write("\n\n")

    f.write("4. Top 30 metastasis-up genes by log2FC\n")
    f.write(top_meta.head(30)[["gene", "log2FC", "score", "padj"]].to_string(index=False))
    f.write("\n\n")

    f.write("5. Top 30 primary-up genes by log2FC\n")
    f.write(top_primary.head(30)[["gene", "log2FC", "score", "padj"]].to_string(index=False))
    f.write("\n\n")

    f.write("6. Files generated\n")
    for x in sorted(os.listdir(FIG5)):
        f.write(x + "\n")

print("DONE FIGURE 5A")
print("Output folder:", FIG5)
print("DE table:", de_out)
print("Report:", report_path)
print("\nCOPY THIS SUMMARY:")
print(open(report_path, encoding="utf-8").read())