import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")

adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")
meta_genes_file = os.path.join(FIG5, "Figure5A_top50_metastasis_up_genes.csv")
primary_genes_file = os.path.join(FIG5, "Figure5A_top50_primary_up_genes.csv")

adata = sc.read_h5ad(adata_path)
adata = adata[adata.obs["label"].isin(["primary", "metastasis"])].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

meta_genes = pd.read_csv(meta_genes_file)["gene"].head(25).tolist()
primary_genes = pd.read_csv(primary_genes_file)["gene"].head(25).tolist()

genes = meta_genes + primary_genes
genes = [g for g in genes if g in adata.var_names]

print("Total genes used:", len(genes))

# FIX: convert categorical columns to string
adata.obs["label_str"] = adata.obs["label"].astype(str)
adata.obs["cell_type_str"] = adata.obs["cell_type"].astype(str)
adata.obs["group"] = adata.obs["label_str"] + "_" + adata.obs["cell_type_str"]

group_order = []
for label in ["primary", "metastasis"]:
    for ct in ["Tumor_epithelial", "Cytotoxic_T_NK", "Plasma_cells", "Macrophage", "Myeloid"]:
        g = label + "_" + ct
        if g in adata.obs["group"].unique():
            group_order.append(g)

expr_rows = []
used_groups = []

for g in group_order:
    sub = adata[adata.obs["group"] == g]
    
    if sub.n_obs < 50:
        continue
    
    mean_expr = sub[:, genes].X.mean(axis=0)
    
    if hasattr(mean_expr, "A1"):
        mean_expr = mean_expr.A1
    else:
        mean_expr = np.asarray(mean_expr).flatten()
    
    expr_rows.append(mean_expr)
    used_groups.append(g)

heatmap_df = pd.DataFrame(expr_rows, index=used_groups, columns=genes)

# Z-score per gene across groups, better for showing gene program differences
heatmap_z = (heatmap_df - heatmap_df.mean(axis=0)) / (heatmap_df.std(axis=0) + 1e-6)

heatmap_df.to_csv(os.path.join(FIG5, "Figure5B_expression_raw.csv"))
heatmap_z.to_csv(os.path.join(FIG5, "Figure5B_expression_zscore.csv"))

plt.figure(figsize=(14, 7))

sns.heatmap(
    heatmap_z,
    cmap="RdBu_r",
    center=0,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "Gene-wise Z-score"}
)

plt.title("Metastasis and primary gene programs across cell types", fontsize=13)
plt.xlabel("Differentially expressed genes")
plt.ylabel("Condition and cell type")

plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=8)

plt.tight_layout()

plt.savefig(
    os.path.join(FIG5, "Figure5B_heatmap.tiff"),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(FIG5, "Figure5B_heatmap.pdf"),
    bbox_inches="tight"
)

plt.close()

report = os.path.join(FIG5, "Figure5B_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:
    f.write("FIGURE 5B REPORT\n\n")
    f.write("Genes used:\n")
    f.write(", ".join(genes) + "\n\n")
    f.write("Groups:\n")
    f.write("\n".join(heatmap_df.index.tolist()) + "\n\n")
    f.write("Shape:\n")
    f.write(str(heatmap_df.shape) + "\n\n")
    f.write("Mean expression table preview:\n")
    f.write(heatmap_df.iloc[:, :10].round(3).to_string())

print("DONE FIGURE 5B")
print(report)