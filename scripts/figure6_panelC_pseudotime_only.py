import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

OUT = r"D:\CRC_META_FULL_SCVI"
FIG6 = os.path.join(OUT, "FIGURE_6_TRAJECTORY")
os.makedirs(FIG6, exist_ok=True)

main_path = os.path.join(
    OUT,
    "FIGURE_5_METASTASIS_MECHANISM",
    "Figure5E_scored_integrated_object.h5ad"
)

adata_full = sc.read_h5ad(main_path)

tumor = adata_full[
    (adata_full.obs["cell_type"].astype(str) == "Tumor_epithelial") &
    (adata_full.obs["label"].astype(str).isin(["primary", "metastasis"]))
].copy()

print("Tumor cells:", tumor.shape)

# Recompute neighbors / pseudotime only
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

sc.tl.diffmap(tumor)

primary_mask = tumor.obs["label"].astype(str) == "primary"
dc1 = tumor.obsm["X_diffmap"][:, 1]

root_candidates = np.where(primary_mask)[0]
root_index = root_candidates[np.argmin(dc1[root_candidates])]

tumor.uns["iroot"] = int(root_index)
sc.tl.dpt(tumor, n_dcs=10)

pt = tumor.obs["dpt_pseudotime"].copy()
if pt[~primary_mask].mean() < pt[primary_mask].mean():
    tumor.obs["dpt_pseudotime"] = 1 - pt

# Save pseudotime object for future
tumor.write(os.path.join(FIG6, "Figure6_tumor_epithelial_trajectory.h5ad"))

df = tumor.obs.copy()
df["label"] = df["label"].astype(str)

plt.figure(figsize=(5.5, 4))

colors = {
    "primary": "#4575B4",
    "metastasis": "#D73027"
}

for lab in ["primary", "metastasis"]:
    vals = df[df["label"] == lab]["dpt_pseudotime"]
    sns.kdeplot(
        vals,
        label=lab,
        fill=True,
        alpha=0.3,
        linewidth=2,
        color=colors[lab]
    )

xmin = df["dpt_pseudotime"].quantile(0.02)
xmax = df["dpt_pseudotime"].quantile(0.98)
plt.xlim(xmin, xmax)

plt.xlabel("Pseudotime")
plt.ylabel("Density")
plt.legend(frameon=False)
sns.despine()
plt.tight_layout()

out_base = os.path.join(FIG6, "C_pseudotime_distribution")
plt.savefig(out_base + ".tiff", dpi=600, bbox_inches="tight")
plt.savefig(out_base + ".pdf", bbox_inches="tight")
plt.close()

print("DONE PANEL C")
print(out_base + ".tiff")
print("Mean pseudotime:")
print(df.groupby("label")["dpt_pseudotime"].mean())