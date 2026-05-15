import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")

de_file = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
de = pd.read_csv(de_file)

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------
de["padj_plot"] = de["padj"].replace(0, 1e-300)
de["neg_log10_padj"] = -np.log10(de["padj_plot"])

# Cap extreme values (visual clarity)
y_cap = 300
de["neg_log10_padj_capped"] = de["neg_log10_padj"].clip(upper=y_cap)

# ------------------------------------------------------------
# Define label genes (clean biological anchors)
# ------------------------------------------------------------
metastasis_genes_to_label = [
    "HSPA1A", "HSPA6", "XCL1", "XCL2", "CCL5", "KLRF1", "IL12RB2", "KLRC1"
]

primary_genes_to_label = [
    "IGFBP6", "CAV1", "CAV2", "COL6A3", "CXCL1", "KRT8", "TM4SF1", "IFI27"
]

label_genes = metastasis_genes_to_label + primary_genes_to_label
label_df = de[de["gene"].isin(label_genes)].copy()

# ------------------------------------------------------------
# Split groups
# ------------------------------------------------------------
ns = de[de["direction"] == "not_significant"]
primary = de[de["direction"] == "primary_up"]
meta = de[de["direction"] == "metastasis_up"]

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8))

# Background (not significant)
ax.scatter(
    ns["log2FC"],
    ns["neg_log10_padj_capped"],
    s=5,
    c="#D0D0D0",
    alpha=0.35,
    linewidths=0
)

# Primary-up (LEFT)
ax.scatter(
    primary["log2FC"],
    primary["neg_log10_padj_capped"],
    s=9,
    c="#4575B4",
    alpha=0.65,
    linewidths=0
)

# Metastasis-up (RIGHT)
ax.scatter(
    meta["log2FC"],
    meta["neg_log10_padj_capped"],
    s=9,
    c="#D73027",
    alpha=0.75,
    linewidths=0
)

# Threshold lines
ax.axvline(-0.25, linestyle="--", linewidth=0.8, color="black", alpha=0.4)
ax.axvline(0.25, linestyle="--", linewidth=0.8, color="black", alpha=0.4)
ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8, color="black", alpha=0.4)

# ------------------------------------------------------------
# Labels (clean placement)
# ------------------------------------------------------------
texts = []

for _, row in label_df.iterrows():
    x = row["log2FC"]
    y = row["neg_log10_padj_capped"]

    # Push labels outward
    if row["gene"] in metastasis_genes_to_label:
        x_text = x + 0.18
        ha = "left"
    else:
        x_text = x - 0.18
        ha = "right"

    texts.append(
        ax.text(
            x_text,
            y,
            row["gene"],
            fontsize=8,
            ha=ha,
            va="center"
        )
    )

adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="black", lw=0.35, alpha=0.5),
    expand_text=(1.3, 1.5),
    force_text=(0.5, 0.7),
    lim=700
)

# ------------------------------------------------------------
# Axes
# ------------------------------------------------------------
ax.set_xlabel("log2 fold change: metastasis vs primary", fontsize=12)
ax.set_ylabel("-log10 adjusted p-value", fontsize=12)

ax.set_title(
    "Primary versus metastasis differential expression",
    fontsize=13,
    pad=25
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_xlim(
    min(de["log2FC"].min() - 0.3, -4.5),
    max(de["log2FC"].max() + 0.3, 3.5)
)

# ------------------------------------------------------------
# TOP HEADER (clean legend replacement)
# ------------------------------------------------------------
fig.text(
    0.02, 0.985,
    f"Primary-up ({len(primary)})",
    fontsize=11,
    color="#4575B4",
    ha="left",
    va="top"
)

fig.text(
    0.50, 0.985,
    f"Not significant ({len(ns)})",
    fontsize=10,
    color="#555555",
    ha="center",
    va="top"
)

fig.text(
    0.98, 0.985,
    f"Metastasis-up ({len(meta)})",
    fontsize=11,
    color="#D73027",
    ha="right",
    va="top"
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
plt.tight_layout(rect=[0, 0, 1, 0.95])

out_tiff = os.path.join(FIG5, "Figure5A_volcano_FINAL.tiff")
out_pdf = os.path.join(FIG5, "Figure5A_volcano_FINAL.pdf")

plt.savefig(out_tiff, dpi=600, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

# Save labeled genes
label_df.to_csv(
    os.path.join(FIG5, "Figure5A_labeled_genes_FINAL.csv"),
    index=False
)

print("DONE FINAL FIGURE 5A")
print(out_tiff)
print(out_pdf)