import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")

de_file = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
de = pd.read_csv(de_file)

# Avoid infinite y values
de["padj_plot"] = de["padj"].replace(0, 1e-300)
de["neg_log10_padj"] = -np.log10(de["padj_plot"])

# Choose cleaner labels:
# top metastasis-up by score + log2FC
meta_label = (
    de[(de["direction"] == "metastasis_up")]
    .sort_values(["score", "log2FC"], ascending=[False, False])
    .head(12)
)

primary_label = (
    de[(de["direction"] == "primary_up")]
    .sort_values(["score", "log2FC"], ascending=[True, True])
    .head(12)
)

label_df = pd.concat([meta_label, primary_label]).drop_duplicates("gene")

# Colors
colors = {
    "metastasis_up": "#D73027",
    "primary_up": "#4575B4",
    "not_significant": "#BDBDBD"
}

fig, ax = plt.subplots(figsize=(9, 8))

for direction, sub in de.groupby("direction"):
    ax.scatter(
        sub["log2FC"],
        sub["neg_log10_padj"],
        s=8 if direction != "not_significant" else 5,
        alpha=0.70 if direction != "not_significant" else 0.35,
        c=colors.get(direction, "#BDBDBD"),
        linewidths=0,
        label=direction.replace("_", " ")
    )

# Threshold lines
ax.axvline(0.25, linestyle="--", linewidth=0.8, color="black", alpha=0.5)
ax.axvline(-0.25, linestyle="--", linewidth=0.8, color="black", alpha=0.5)
ax.axhline(-np.log10(0.05), linestyle="--", linewidth=0.8, color="black", alpha=0.5)

# Labels with repulsion
texts = []
for _, row in label_df.iterrows():
    texts.append(
        ax.text(
            row["log2FC"],
            row["neg_log10_padj"],
            row["gene"],
            fontsize=8,
            ha="center",
            va="center"
        )
    )

adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="black", lw=0.4, alpha=0.6),
    expand_text=(1.2, 1.4),
    expand_points=(1.2, 1.4),
    force_text=(0.4, 0.6),
    force_points=(0.2, 0.4),
    lim=500
)

ax.set_xlabel("log2 fold change: metastasis vs primary", fontsize=12)
ax.set_ylabel("-log10 adjusted p-value", fontsize=12)
ax.set_title("Primary versus metastasis differential expression", fontsize=13)

ax.legend(
    frameon=False,
    loc="upper right",
    fontsize=9
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_tiff = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_volcano_CLEAN_labels.tiff")
out_pdf = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_volcano_CLEAN_labels.pdf")

plt.savefig(out_tiff, dpi=600, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

# Save labeled genes
label_df.to_csv(
    os.path.join(FIG5, "Figure5A_labeled_genes_clean_volcano.csv"),
    index=False
)

print("DONE CLEAN VOLCANO")
print(out_tiff)
print(out_pdf)