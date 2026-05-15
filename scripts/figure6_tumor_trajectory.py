import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG6 = os.path.join(OUT, "FIGURE_6_TRAJECTORY")
os.makedirs(FIG6, exist_ok=True)

adata_path = os.path.join(
    OUT,
    "FIGURE_5_METASTASIS_MECHANISM",
    "Figure5E_scored_integrated_object.h5ad"
)

if not os.path.exists(adata_path):
    adata_path = os.path.join(OUT, "05_CRC_CORE_scVI_integrated_celltype_annotated.h5ad")

adata = sc.read_h5ad(adata_path)

# ------------------------------------------------------------
# Subset tumor epithelial
# ------------------------------------------------------------
tumor = adata[
    (adata.obs["cell_type"].astype(str) == "Tumor_epithelial") &
    (adata.obs["label"].astype(str).isin(["primary", "metastasis"]))
].copy()

print("Tumor object:", tumor)

# ------------------------------------------------------------
# Neighbors / embedding
# ------------------------------------------------------------
if "X_scVI" in tumor.obsm:
    sc.pp.neighbors(tumor, use_rep="X_scVI", n_neighbors=25)
else:
    sc.pp.normalize_total(tumor, target_sum=1e4)
    sc.pp.log1p(tumor)
    sc.pp.highly_variable_genes(tumor, n_top_genes=2000)
    tumor = tumor[:, tumor.var["highly_variable"]].copy()
    sc.pp.scale(tumor, max_value=10)
    sc.tl.pca(tumor, n_comps=30)
    sc.pp.neighbors(tumor, use_rep="X_pca", n_neighbors=25)

sc.tl.umap(tumor, min_dist=0.35)
sc.tl.diffmap(tumor)

# ------------------------------------------------------------
# Root = primary cells
# ------------------------------------------------------------
primary_mask = tumor.obs["label"].astype(str) == "primary"
dc1 = tumor.obsm["X_diffmap"][:, 1]

root_candidates = np.where(primary_mask)[0]
root_index = root_candidates[np.argmin(dc1[root_candidates])]

tumor.uns["iroot"] = int(root_index)

sc.tl.dpt(tumor, n_dcs=10)

# Orient pseudotime
pt = tumor.obs["dpt_pseudotime"].copy()
mean_primary = pt[primary_mask].mean()
mean_meta = pt[~primary_mask].mean()

if mean_meta < mean_primary:
    tumor.obs["dpt_pseudotime"] = 1 - tumor.obs["dpt_pseudotime"]

# ------------------------------------------------------------
# Save helper
# ------------------------------------------------------------
def save_umap(color, filename, title):
    sc.pl.umap(
        tumor,
        color=color,
        frameon=False,
        size=4,
        show=False,
        title=title
    )
    plt.savefig(os.path.join(FIG6, filename + ".tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG6, filename + ".pdf"), bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# Basic plots
# ------------------------------------------------------------
save_umap("label", "Figure6A_label", "Primary vs Metastasis")
save_umap("dpt_pseudotime", "Figure6B_pseudotime", "Pseudotime")

# ------------------------------------------------------------
# Pathway scores
# ------------------------------------------------------------
scores = [
    "Metastasis_immune_score",
    "Metastasis_stress_score",
    "Primary_structural_score",
    "Primary_ECM_EMT_score"
]

for score in scores:
    if score in tumor.obs.columns:
        title = score.replace("_", " ")
        filename = "Figure6_score_" + score
        save_umap(score, filename, title)

# ------------------------------------------------------------
# Histogram
# ------------------------------------------------------------
df = tumor.obs.copy()
df["label"] = df["label"].astype(str)

plt.figure(figsize=(7,5))

for lab, color in [("primary","#4575B4"), ("metastasis","#D73027")]:
    vals = df[df["label"] == lab]["dpt_pseudotime"]
    plt.hist(vals, bins=60, density=True, alpha=0.5, label=lab, color=color)

plt.xlabel("Pseudotime")
plt.ylabel("Density")
plt.legend(frameon=False)
plt.title("Pseudotime distribution")
plt.tight_layout()

plt.savefig(os.path.join(FIG6, "Figure6C_hist.tiff"), dpi=600)
plt.savefig(os.path.join(FIG6, "Figure6C_hist.pdf"))
plt.close()

# ------------------------------------------------------------
# Report
# ------------------------------------------------------------
report = os.path.join(FIG6, "Figure6_TRAJECTORY_REPORT.txt")

with open(report, "w") as f:
    f.write("FIGURE 6 TRAJECTORY REPORT\n\n")
    f.write("Cells:\n")
    f.write(str(tumor.shape) + "\n\n")
    f.write("Mean pseudotime by label:\n")
    f.write(df.groupby("label")["dpt_pseudotime"].mean().to_string())

print("DONE FIGURE 6")
print(report)