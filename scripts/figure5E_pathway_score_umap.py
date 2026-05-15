import os
import scanpy as sc
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
os.makedirs(FIG5, exist_ok=True)

adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")
adata = sc.read_h5ad(adata_path)

# Normalize for gene scoring
adata_expr = adata.copy()
sc.pp.normalize_total(adata_expr, target_sum=1e4)
sc.pp.log1p(adata_expr)

# Gene programs based on your Figure 5A-D results
gene_sets = {
    "Metastasis_immune_score": [
        "XCL1", "XCL2", "CCL5", "IFNG", "TNF", "GZMA", "GZMK",
        "GNLY", "PRF1", "KLRF1", "KLRC1", "KLRC2", "CD160",
        "IL12RB2", "TNFSF13", "TNFSF14"
    ],
    "Metastasis_stress_score": [
        "HSPA1A", "HSPA1B", "HSPA6", "HSPH1", "HSPD1",
        "HSPA8", "DNAJB1", "DNAJA4", "HSPE1", "BAG3"
    ],
    "Primary_structural_score": [
        "IGFBP6", "CAV1", "CAV2", "COL6A3", "COL6A1",
        "COL5A1", "F3", "MGP", "TM4SF1", "KRT8", "KRT18",
        "ITGB8", "CXCL1", "RARRES2"
    ],
    "Primary_ECM_EMT_score": [
        "COL6A3", "COL6A1", "COL5A1", "CAV1", "CAV2",
        "MGP", "ITGB8", "PLAU", "TGM2", "JUP", "TJP1",
        "LAMB2", "LAMB3"
    ]
}

# Score genes
used = {}

for score_name, genes in gene_sets.items():
    valid = [g for g in genes if g in adata_expr.var_names]
    used[score_name] = valid
    print(score_name, "valid genes:", valid)

    if len(valid) >= 3:
        sc.tl.score_genes(
            adata_expr,
            gene_list=valid,
            score_name=score_name,
            use_raw=False
        )
    else:
        print("SKIPPED", score_name, "not enough genes")

# Transfer scores to original adata
for score_name in used.keys():
    if score_name in adata_expr.obs.columns:
        adata.obs[score_name] = adata_expr.obs[score_name]

# Save updated object
adata.write(os.path.join(FIG5, "Figure5E_scored_integrated_object.h5ad"))

# Save score table
score_cols = [s for s in gene_sets.keys() if s in adata.obs.columns]
adata.obs[
    ["dataset", "label", "site", "cell_type", "cluster"] + score_cols
].to_csv(os.path.join(FIG5, "Figure5E_pathway_scores_metadata.csv"))

# Plot each score separately
for score in score_cols:
    sc.pl.umap(
        adata,
        color=score,
        cmap="viridis",
        frameon=False,
        size=2,
        show=False,
        title=score.replace("_", " ")
    )

    plt.savefig(
        os.path.join(FIG5, f"Figure5E_{score}_UMAP.tiff"),
        dpi=600,
        bbox_inches="tight"
    )

    plt.savefig(
        os.path.join(FIG5, f"Figure5E_{score}_UMAP.pdf"),
        bbox_inches="tight"
    )

    plt.close()

# Also make label/celltype reference UMAPs
for color in ["label", "cell_type"]:
    sc.pl.umap(
        adata,
        color=color,
        frameon=False,
        size=2,
        show=False,
        title=color
    )

    plt.savefig(
        os.path.join(FIG5, f"Figure5E_reference_{color}_UMAP.tiff"),
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()

# Summary by label and cell type
import pandas as pd

summary_label = adata.obs.groupby("label", observed=True)[score_cols].mean()
summary_celltype = adata.obs.groupby("cell_type", observed=True)[score_cols].mean()
summary_label_celltype = adata.obs.groupby(["label", "cell_type"], observed=True)[score_cols].mean()

summary_label.to_csv(os.path.join(FIG5, "Figure5E_score_summary_by_label.csv"))
summary_celltype.to_csv(os.path.join(FIG5, "Figure5E_score_summary_by_celltype.csv"))
summary_label_celltype.to_csv(os.path.join(FIG5, "Figure5E_score_summary_by_label_celltype.csv"))

report_path = os.path.join(FIG5, "Figure5E_PATHWAY_SCORE_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 5E PATHWAY SCORE UMAP REPORT\n\n")

    f.write("1. Gene sets used\n")
    for k, v in used.items():
        f.write(k + ": " + ", ".join(v) + "\n")
    f.write("\n")

    f.write("2. Score summary by label\n")
    f.write(summary_label.to_string())
    f.write("\n\n")

    f.write("3. Score summary by cell type\n")
    f.write(summary_celltype.to_string())
    f.write("\n\n")

    f.write("4. Score summary by label and cell type\n")
    f.write(summary_label_celltype.to_string())
    f.write("\n\n")

    f.write("5. Files generated\n")
    for x in sorted(os.listdir(FIG5)):
        if "Figure5E" in x:
            f.write(x + "\n")

print("DONE FIGURE 5E")
print(report_path)
print(open(report_path, encoding="utf-8").read())