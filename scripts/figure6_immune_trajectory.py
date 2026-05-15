import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG6 = os.path.join(OUT, "FIGURE_6_TRAJECTORY_IMMUNE")
os.makedirs(FIG6, exist_ok=True)

adata_path = os.path.join(
    OUT,
    "FIGURE_5_METASTASIS_MECHANISM",
    "Figure5E_scored_integrated_object.h5ad"
)

if not os.path.exists(adata_path):
    adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")

adata = sc.read_h5ad(adata_path)

immune_types = ["Cytotoxic_T_NK", "Macrophage", "Myeloid", "Plasma_cells"]

immune = adata[
    (adata.obs["cell_type"].astype(str).isin(immune_types)) &
    (adata.obs["label"].astype(str).isin(["primary", "metastasis"]))
].copy()

print("Immune object:", immune)

if "X_scVI" in immune.obsm:
    sc.pp.neighbors(immune, use_rep="X_scVI", n_neighbors=30)
else:
    sc.pp.normalize_total(immune, target_sum=1e4)
    sc.pp.log1p(immune)
    sc.pp.highly_variable_genes(immune, n_top_genes=2000)
    immune = immune[:, immune.var["highly_variable"]].copy()
    sc.pp.scale(immune)
    sc.tl.pca(immune)
    sc.pp.neighbors(immune)

sc.tl.umap(immune, min_dist=0.3)
sc.tl.diffmap(immune)

primary_mask = immune.obs["label"].astype(str) == "primary"
dc1 = immune.obsm["X_diffmap"][:, 1]

root_candidates = np.where(primary_mask)[0]
root_index = root_candidates[np.argmin(dc1[root_candidates])]

immune.uns["iroot"] = int(root_index)
sc.tl.dpt(immune)

pt = immune.obs["dpt_pseudotime"]
if pt[~primary_mask].mean() < pt[primary_mask].mean():
    immune.obs["dpt_pseudotime"] = 1 - pt

immune.write(os.path.join(FIG6, "Figure6_immune_trajectory.h5ad"))

def save(color, name):
    sc.pl.umap(
        immune,
        color=color,
        frameon=False,
        size=3,
        show=False,
        title="",
        legend_loc="right margin",
        legend_fontsize=10
    )

    plt.subplots_adjust(right=0.75)

    plt.savefig(
        os.path.join(FIG6, name + ".tiff"),
        dpi=600,
        bbox_inches="tight"
    )

    plt.savefig(
        os.path.join(FIG6, name + ".pdf"),
        bbox_inches="tight"
    )

    plt.close()

save("label", "A_label")
save("dpt_pseudotime", "B_pseudotime")
save("cell_type", "C_celltype")

for score in [
    "Metastasis_immune_score",
    "Metastasis_stress_score",
    "Primary_structural_score"
]:
    if score in immune.obs.columns:
        save(score, "D_" + score)

report = os.path.join(FIG6, "Figure6_IMMUNE_REPORT.txt")
df = immune.obs.copy()

with open(report, "w", encoding="utf-8") as f:
    f.write("IMMUNE TRAJECTORY REPORT\n\n")
    f.write("Cells:\n")
    f.write(str(immune.shape) + "\n\n")
    f.write("Mean pseudotime by label:\n")
    f.write(df.groupby("label", observed=True)["dpt_pseudotime"].mean().to_string())
    f.write("\n\n")
    f.write("Mean pseudotime by cell type:\n")
    f.write(df.groupby("cell_type", observed=True)["dpt_pseudotime"].mean().to_string())

print("DONE IMMUNE TRAJECTORY")
print(report)
print(open(report, encoding="utf-8").read())