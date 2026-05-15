import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

# -----------------------------
# Panel A: signature overview
# -----------------------------
meta_file = os.path.join(FIG7, "7A_metastasis_signature_top300.csv")
primary_file = os.path.join(FIG7, "7A_primary_signature_top300.csv")
cand_file = os.path.join(FIG7, "7B_top50_candidate_drug_reversal.csv")

meta = pd.read_csv(meta_file)
primary = pd.read_csv(primary_file)

# simple top signature display
meta_genes = meta.iloc[:25, 0].astype(str).tolist()
primary_genes = primary.iloc[:25, 0].astype(str).tolist()

fig, ax = plt.subplots(figsize=(8, 7))
ax.axis("off")

ax.text(0.02, 0.98, "Metastasis-up signature", fontsize=14, fontweight="bold", va="top", color="#D73027")
ax.text(0.52, 0.98, "Primary-up signature", fontsize=14, fontweight="bold", va="top", color="#4575B4")

for i, g in enumerate(meta_genes):
    ax.text(0.02, 0.92 - i * 0.035, g, fontsize=10, va="top", color="#D73027")

for i, g in enumerate(primary_genes):
    ax.text(0.52, 0.92 - i * 0.035, g, fontsize=10, va="top", color="#4575B4")

plt.tight_layout()
plt.savefig(os.path.join(FIG7, "7A_signature_gene_overview.tiff"), dpi=600, bbox_inches="tight")
plt.savefig(os.path.join(FIG7, "7A_signature_gene_overview.pdf"), bbox_inches="tight")
plt.close()

# -----------------------------
# Panel C: clean bubble plot
# -----------------------------
candidates = pd.read_csv(cand_file).head(30).copy()

candidates["minus_log10_adj_p"] = -np.log10(candidates["Adjusted P-value"].clip(lower=1e-300))
candidates["rank"] = range(1, len(candidates) + 1)

# extract overlap numerator
candidates["overlap_n"] = (
    candidates["Overlap"]
    .astype(str)
    .str.extract(r"(\d+)")[0]
    .fillna(5)
    .astype(float)
)

# save ranked map for figure legend
map_cols = ["rank", "Term", "library", "Adjusted P-value", "Combined Score", "Overlap", "Genes", "reversal_class", "rank_score"]
map_cols = [c for c in map_cols if c in candidates.columns]
candidates[map_cols].to_csv(os.path.join(FIG7, "7C_bubble_rank_to_drug_name_table.csv"), index=False)

plt.figure(figsize=(8, 6))

plt.scatter(
    candidates["Combined Score"],
    candidates["minus_log10_adj_p"],
    s=np.clip(candidates["overlap_n"] * 18, 60, 450),
    alpha=0.65,
    c="#6A3D9A",
    edgecolor="black",
    linewidth=0.3
)

# Use rank numbers instead of full drug names
for _, row in candidates.iterrows():
    plt.text(
        row["Combined Score"],
        row["minus_log10_adj_p"],
        str(int(row["rank"])),
        fontsize=7,
        ha="center",
        va="center",
        color="white",
        fontweight="bold"
    )

plt.xlabel("Combined score")
plt.ylabel("-log10 adjusted p-value")
plt.title("")
plt.tight_layout()

plt.savefig(os.path.join(FIG7, "7C_drug_significance_bubble_CLEAN.tiff"), dpi=600, bbox_inches="tight")
plt.savefig(os.path.join(FIG7, "7C_drug_significance_bubble_CLEAN.pdf"), bbox_inches="tight")
plt.close()

# -----------------------------
# Panel C side table for top 12
# -----------------------------
top12 = candidates.head(12).copy()
top12["short_term"] = top12["Term"].astype(str).apply(lambda x: "\n".join(textwrap.wrap(x, width=34)))

fig, ax = plt.subplots(figsize=(9, 6))
ax.axis("off")

y = 0.98
ax.text(0.02, y, "Rank", fontsize=11, fontweight="bold", va="top")
ax.text(0.14, y, "Candidate perturbation / drug signature", fontsize=11, fontweight="bold", va="top")
ax.text(0.78, y, "Overlap", fontsize=11, fontweight="bold", va="top")

y -= 0.07
for _, row in top12.iterrows():
    ax.text(0.02, y, str(int(row["rank"])), fontsize=9, va="top")
    ax.text(0.14, y, row["short_term"], fontsize=8, va="top")
    ax.text(0.78, y, str(row["Overlap"]), fontsize=8, va="top")
    y -= 0.075

plt.tight_layout()
plt.savefig(os.path.join(FIG7, "7C_top12_candidate_name_table.tiff"), dpi=600, bbox_inches="tight")
plt.savefig(os.path.join(FIG7, "7C_top12_candidate_name_table.pdf"), bbox_inches="tight")
plt.close()

print("DONE fixed panels")
print("Created:")
print("7A_signature_gene_overview.tiff")
print("7C_drug_significance_bubble_CLEAN.tiff")
print("7C_top12_candidate_name_table.tiff")
print("7C_bubble_rank_to_drug_name_table.csv")