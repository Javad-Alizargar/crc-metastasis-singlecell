import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_REBUILD")

file = os.path.join(FIG7, "7A_metastasis_signature_all_genes_annotated.csv")
df = pd.read_csv(file)

# focus on top 120 (important)
df = df.sort_values("rank_log2FC").head(120).copy()

# colors per module
colors = {
    "Stress_heatshock": "#E41A1C",
    "Cytokine_inflammatory": "#FF7F00",
    "Cytotoxic_NK_T": "#377EB8",
    "Signaling_regulatory": "#984EA3",
    "Metabolic_lipid": "#4DAF4A",
    "Antigen_presentation_myeloid": "#A65628",
    "Other": "#BDBDBD"
}

df["color"] = df["module"].map(colors).fillna("#BDBDBD")

# size scaling
df["size"] = np.clip(df["minus_log10_padj"], 1, 300)

# normalize for plotting
df["size"] = df["size"] / df["size"].max() * 120

plt.figure(figsize=(10,6))

# scatter
plt.scatter(
    df["rank_log2FC"],
    df["log2FC"],
    c=df["color"],
    s=df["size"],
    alpha=0.8,
    edgecolor="black",
    linewidth=0.2
)

# vertical lollipops
for _, row in df.head(30).iterrows():
    plt.plot(
        [row["rank_log2FC"], row["rank_log2FC"]],
        [0, row["log2FC"]],
        color=row["color"],
        alpha=0.3,
        linewidth=1
    )

# label only key genes (avoid clutter)
key = df[
    (df["rank_log2FC"] <= 25) |
    (df["module"].isin(["Stress_heatshock", "Cytokine_inflammatory"]))
].head(35)

for _, row in key.iterrows():
    plt.text(
        row["rank_log2FC"],
        row["log2FC"] + 0.05,
        row["gene"],
        fontsize=7,
        rotation=45
    )

plt.xlabel("Gene rank (metastasis-up)")
plt.ylabel("log2 fold-change")

plt.xlim(0, 120)
plt.ylim(0, df["log2FC"].max() + 0.5)

plt.title("")

plt.tight_layout()

plt.savefig(os.path.join(FIG7, "7A_signature_landscape.tiff"), dpi=600)
plt.savefig(os.path.join(FIG7, "7A_signature_landscape.pdf"))

plt.close()

print("DONE 7A PLOT")