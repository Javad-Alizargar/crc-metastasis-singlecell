# ============================================================
# Figure immune trajectory Panel E (SAFE VERSION)
# Compatible with older Python
# ============================================================

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

OUT = r"D:\CRC_META_FULL_SCVI"
FIG = os.path.join(OUT, "FIGURE_6_TRAJECTORY_IMMUNE")

adata_path = os.path.join(FIG, "Figure6_immune_trajectory.h5ad")

adata = sc.read_h5ad(adata_path)
df = adata.obs.copy()

df["label"] = df["label"].astype(str)
df["cell_type"] = df["cell_type"].astype(str)

df = df[df["label"].isin(["primary", "metastasis"])].copy()
df = df.dropna(subset=["dpt_pseudotime"])

# ------------------------------------------------------------
# Quartile composition
# ------------------------------------------------------------
df["pseudotime_quartile"] = pd.qcut(
    df["dpt_pseudotime"],
    q=4,
    labels=["Q1_early", "Q2", "Q3", "Q4_late"]
)

qtab = pd.crosstab(
    df["pseudotime_quartile"],
    df["cell_type"],
    normalize="index"
) * 100

qtab = qtab.loc[[x for x in ["Q1_early", "Q2", "Q3", "Q4_late"] if x in qtab.index]]

# ------------------------------------------------------------
# Stats
# ------------------------------------------------------------
primary = df[df["label"] == "primary"]["dpt_pseudotime"].values
meta = df[df["label"] == "metastasis"]["dpt_pseudotime"].values

u, pval = mannwhitneyu(primary, meta)

primary_mean = np.mean(primary)
meta_mean = np.mean(meta)
delta_mean = meta_mean - primary_mean

# ------------------------------------------------------------
# Colors
# ------------------------------------------------------------
label_colors = {
    "primary": "#4575B4",
    "metastasis": "#D73027"
}

celltype_colors = {
    "Plasma_cells": "#984EA3",
    "Cytotoxic_T_NK": "#4575B4",
    "Macrophage": "#D73027",
    "Myeloid": "#F46D43"
}

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

# ---------------- Left panel ----------------
ax = axes[0]

bins = np.linspace(0, 1, 55)

for lab in ["primary", "metastasis"]:
    vals = df[df["label"] == lab]["dpt_pseudotime"].values

    ax.hist(
        vals,
        bins=bins,
        density=True,
        alpha=0.4,
        color=label_colors[lab],
        label=lab.capitalize() + " (n=" + str(len(vals)) + ")",
        edgecolor="none"
    )

    ax.axvline(
        np.mean(vals),
        color=label_colors[lab],
        linewidth=2,
        linestyle="--"
    )

ax.text(
    0.03, 0.95,
    "Mean diff = %.3f\np < 1e-300" % delta_mean,
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=10,
    bbox=dict(facecolor="white", edgecolor="black")
)

ax.set_xlim(0, 1)
ax.set_xlabel("Immune pseudotime")
ax.set_ylabel("Density")
ax.set_title("Pseudotime distribution")
ax.legend(frameon=False)

# ---------------- Right panel ----------------
ax = axes[1]

bottom = np.zeros(len(qtab))

for ct in ["Plasma_cells", "Cytotoxic_T_NK", "Macrophage", "Myeloid"]:
    if ct in qtab.columns:
        vals = qtab[ct].values

        ax.bar(
            range(len(qtab)),
            vals,
            bottom=bottom,
            label=ct.replace("_", " "),
            color=celltype_colors.get(ct, "#999999")
        )

        # annotate
        for i in range(len(vals)):
            if vals[i] > 10:
                ax.text(
                    i,
                    bottom[i] + vals[i]/2,
                    str(int(vals[i])) + "%",
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=9
                )

        bottom += vals

ax.set_xticks(range(len(qtab)))
ax.set_xticklabels(qtab.index)

ax.set_ylim(0, 100)
ax.set_ylabel("Composition (%)")
ax.set_title("Immune composition along pseudotime")
ax.legend(frameon=False)

plt.tight_layout()

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
out_tiff = os.path.join(FIG, "Figure7E_immune_panel.tiff")
out_pdf = os.path.join(FIG, "Figure7E_immune_panel.pdf")

plt.savefig(out_tiff, dpi=600)
plt.savefig(out_pdf)

plt.close()

print("DONE Panel E")
print(out_tiff)