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

adata = adata[adata.obs["label"].isin(["primary", "metastasis"])].copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

cell_types = [
    "Tumor_epithelial",
    "Cytotoxic_T_NK",
    "Plasma_cells",
    "Macrophage",
    "Myeloid"
]

all_results = []
summary_rows = []

for ct in cell_types:
    print("\nRunning DE for:", ct)

    sub = adata[adata.obs["cell_type"].astype(str) == ct].copy()

    label_counts = sub.obs["label"].value_counts().to_dict()
    n_primary = label_counts.get("primary", 0)
    n_meta = label_counts.get("metastasis", 0)

    print("primary:", n_primary, "metastasis:", n_meta)

    if n_primary < 100 or n_meta < 100:
        summary_rows.append({
            "cell_type": ct,
            "primary_cells": n_primary,
            "metastasis_cells": n_meta,
            "metastasis_up": None,
            "primary_up": None,
            "not_significant": None,
            "status": "skipped_low_cells"
        })
        continue

    sc.tl.rank_genes_groups(
        sub,
        groupby="label",
        groups=["metastasis"],
        reference="primary",
        method="wilcoxon",
        n_genes=sub.n_vars
    )

    de = sc.get.rank_genes_groups_df(sub, group="metastasis")

    de = de.rename(columns={
        "names": "gene",
        "scores": "score",
        "logfoldchanges": "log2FC",
        "pvals": "pval",
        "pvals_adj": "padj"
    })

    de["cell_type"] = ct
    de["primary_cells"] = n_primary
    de["metastasis_cells"] = n_meta
    de["neg_log10_padj"] = -np.log10(de["padj"].replace(0, 1e-300))

    de["direction"] = "not_significant"
    de.loc[(de["padj"] < 0.05) & (de["log2FC"] > 0.25), "direction"] = "metastasis_up"
    de.loc[(de["padj"] < 0.05) & (de["log2FC"] < -0.25), "direction"] = "primary_up"

    de.to_csv(
        os.path.join(FIG5, f"Figure5C_DE_{ct}.csv"),
        index=False
    )

    counts = de["direction"].value_counts().to_dict()

    summary_rows.append({
        "cell_type": ct,
        "primary_cells": n_primary,
        "metastasis_cells": n_meta,
        "metastasis_up": counts.get("metastasis_up", 0),
        "primary_up": counts.get("primary_up", 0),
        "not_significant": counts.get("not_significant", 0),
        "status": "ok"
    })

    top_meta = de[de["direction"] == "metastasis_up"].sort_values("log2FC", ascending=False).head(20)
    top_primary = de[de["direction"] == "primary_up"].sort_values("log2FC", ascending=True).head(20)

    top_meta.to_csv(os.path.join(FIG5, f"Figure5C_top20_metastasis_up_{ct}.csv"), index=False)
    top_primary.to_csv(os.path.join(FIG5, f"Figure5C_top20_primary_up_{ct}.csv"), index=False)

    all_results.append(de)

summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(FIG5, "Figure5C_celltype_DE_summary.csv"), index=False)

if len(all_results) > 0:
    all_de = pd.concat(all_results, ignore_index=True)
    all_de.to_csv(os.path.join(FIG5, "Figure5C_all_celltype_DE_results.csv"), index=False)

# ------------------------------------------------------------
# Panel 5C plot: DE gene counts per cell type
# ------------------------------------------------------------
plot_df = summary[summary["status"] == "ok"].copy()

fig, ax = plt.subplots(figsize=(9, 5))

x = np.arange(len(plot_df))
width = 0.35

ax.bar(
    x - width / 2,
    plot_df["primary_up"],
    width,
    label="Primary-up",
    color="#4575B4"
)

ax.bar(
    x + width / 2,
    plot_df["metastasis_up"],
    width,
    label="Metastasis-up",
    color="#D73027"
)

ax.set_xticks(x)
ax.set_xticklabels(plot_df["cell_type"], rotation=35, ha="right")
ax.set_ylabel("Number of DE genes")
ax.set_title("Cell-type-specific differential expression")
ax.legend(frameon=False)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    os.path.join(FIG5, "Figure5C_celltype_DE_gene_counts.tiff"),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(FIG5, "Figure5C_celltype_DE_gene_counts.pdf"),
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
report_path = os.path.join(FIG5, "Figure5C_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 5C CELL-TYPE-SPECIFIC DE REPORT\n\n")

    f.write("1. Summary table\n")
    f.write(summary.to_string(index=False))
    f.write("\n\n")

    if len(all_results) > 0:
        f.write("2. Top metastasis-up genes per cell type\n")
        for ct in cell_types:
            path = os.path.join(FIG5, f"Figure5C_top20_metastasis_up_{ct}.csv")
            if os.path.exists(path):
                tmp = pd.read_csv(path)
                f.write("\n" + ct + "\n")
                f.write(tmp[["gene", "log2FC", "score", "padj"]].head(10).to_string(index=False))
                f.write("\n")

        f.write("\n\n3. Top primary-up genes per cell type\n")
        for ct in cell_types:
            path = os.path.join(FIG5, f"Figure5C_top20_primary_up_{ct}.csv")
            if os.path.exists(path):
                tmp = pd.read_csv(path)
                f.write("\n" + ct + "\n")
                f.write(tmp[["gene", "log2FC", "score", "padj"]].head(10).to_string(index=False))
                f.write("\n")

print("DONE FIGURE 5C")
print(report_path)
print(open(report_path, encoding="utf-8").read())