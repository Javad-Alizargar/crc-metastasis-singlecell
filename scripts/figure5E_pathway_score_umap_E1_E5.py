import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
os.makedirs(FIG5, exist_ok=True)

adata_path = os.path.join(FIG5, "Figure5E_scored_integrated_object.h5ad")

if not os.path.exists(adata_path):
    adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")

adata = sc.read_h5ad(adata_path)

# If scores are missing, calculate them
score_sets = {
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

missing_scores = [s for s in score_sets if s not in adata.obs.columns]

if missing_scores:
    adata_expr = adata.copy()
    sc.pp.normalize_total(adata_expr, target_sum=1e4)
    sc.pp.log1p(adata_expr)

    for score_name, genes in score_sets.items():
        valid = [g for g in genes if g in adata_expr.var_names]
        if len(valid) >= 3:
            sc.tl.score_genes(
                adata_expr,
                gene_list=valid,
                score_name=score_name,
                use_raw=False
            )
            adata.obs[score_name] = adata_expr.obs[score_name]

    adata.write(os.path.join(FIG5, "Figure5E_scored_integrated_object.h5ad"))

panel_map = [
    ("E1", "Metastasis_immune_score", "metastasis_immune_score_UMAP"),
    ("E2", "Metastasis_stress_score", "metastasis_stress_score_UMAP"),
    ("E3", "Primary_structural_score", "primary_structural_score_UMAP"),
    ("E4", "Primary_ECM_EMT_score", "primary_ECM_EMT_score_UMAP"),
    ("E5", "cell_type", "reference_cell_type_UMAP"),
]

for panel, color, clean_name in panel_map:
    sc.pl.umap(
        adata,
        color=color,
        cmap="viridis" if color.endswith("_score") else None,
        frameon=False,
        size=2,
        show=False,
        title=""
    )

    out_base = os.path.join(FIG5, f"{panel}_{clean_name}")

    plt.savefig(out_base + ".tiff", dpi=600, bbox_inches="tight")
    plt.savefig(out_base + ".pdf", bbox_inches="tight")
    plt.close()

score_cols = [x for x in score_sets if x in adata.obs.columns]

summary_label = adata.obs.groupby("label", observed=True)[score_cols].mean()
summary_celltype = adata.obs.groupby("cell_type", observed=True)[score_cols].mean()
summary_label_celltype = adata.obs.groupby(["label", "cell_type"], observed=True)[score_cols].mean()

summary_label.to_csv(os.path.join(FIG5, "E1_E4_score_summary_by_label.csv"))
summary_celltype.to_csv(os.path.join(FIG5, "E1_E4_score_summary_by_celltype.csv"))
summary_label_celltype.to_csv(os.path.join(FIG5, "E1_E4_score_summary_by_label_celltype.csv"))

report_path = os.path.join(FIG5, "Figure5E_E1_E5_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 5E E1-E5 PATHWAY SCORE REPORT\n\n")
    f.write("Panels generated:\n")
    for panel, color, clean_name in panel_map:
        f.write(f"{panel}: {clean_name}\n")
    f.write("\nScore summary by label:\n")
    f.write(summary_label.to_string())
    f.write("\n\nScore summary by cell type:\n")
    f.write(summary_celltype.to_string())
    f.write("\n\nScore summary by label and cell type:\n")
    f.write(summary_label_celltype.to_string())

print("DONE Figure 5E E1-E5")
print(report_path)