import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

file = os.path.join(FIG7, "7B_candidates_ranked_full.csv")
df = pd.read_csv(file)

# Use top 25 for readability
plot = df.head(25).copy()

# Short readable labels
def short_label(x, width=38):
    x = str(x)
    x = x.replace("LINCS_L1000_Chem_Pert_", "")
    x = x.replace("LJP006 ", "")
    x = x.replace("LJP005 ", "")
    x = x.replace("24H-", "")
    x = x.replace("3H-", "")
    return "\n".join(textwrap.wrap(x, width=width))

plot["plot_label"] = plot["clean_term"].apply(short_label)

# Reverse order so best is at top
plot = plot.iloc[::-1].reset_index(drop=True)

# Effect-class colors
effect_colors = {
    "candidate_reversal": "#1A9850",
    "potential_risk": "#D73027",
    "association": "#6A3D9A"
}
plot["effect_color"] = plot["effect_class"].map(effect_colors).fillna("#999999")

# Module columns
module_cols = [
    "module_Stress_heatshock",
    "module_Signaling_regulatory",
    "module_Other",
    "module_Metabolic_lipid",
    "module_Antigen_presentation_myeloid",
    "module_Cytotoxic_NK_T",
    "module_Cytokine_inflammatory"
]
module_cols = [c for c in module_cols if c in plot.columns]

module_names = [c.replace("module_", "") for c in module_cols]

module_colors = {
    "Stress_heatshock": "#E41A1C",
    "Signaling_regulatory": "#984EA3",
    "Other": "#BDBDBD",
    "Metabolic_lipid": "#4DAF4A",
    "Antigen_presentation_myeloid": "#A65628",
    "Cytotoxic_NK_T": "#377EB8",
    "Cytokine_inflammatory": "#FF7F00"
}

# Figure layout
fig = plt.figure(figsize=(15, 10))
gs = GridSpec(1, 2, width_ratios=[2.2, 1.2], wspace=0.08)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

y = np.arange(len(plot))

# ----------------------------------------
# Left panel: ranked perturbation score
# ----------------------------------------
ax1.hlines(
    y,
    0,
    plot["final_rank_score"],
    color="#CCCCCC",
    linewidth=1.2,
    zorder=1
)

ax1.scatter(
    plot["final_rank_score"],
    y,
    s=np.clip(plot["n_overlap"] * 18, 60, 600),
    c=plot["effect_color"],
    edgecolor="black",
    linewidth=0.4,
    alpha=0.85,
    zorder=2
)

ax1.set_yticks(y)
ax1.set_yticklabels(plot["plot_label"], fontsize=8)
ax1.set_xlabel("Perturbation rank score")
ax1.set_title("Ranked perturbation landscape", fontsize=13)

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Add small text labels for effect class
for i, row in plot.iterrows():
    ax1.text(
        row["final_rank_score"] + plot["final_rank_score"].max() * 0.015,
        i,
        f"{row['direction']} | n={int(row['n_overlap'])}",
        va="center",
        fontsize=7,
        color="#444444"
    )

# ----------------------------------------
# Right panel: module burden stacked bars
# ----------------------------------------
left = np.zeros(len(plot))

for col, name in zip(module_cols, module_names):
    vals = plot[col].fillna(0).values
    ax2.barh(
        y,
        vals,
        left=left,
        color=module_colors.get(name, "#999999"),
        edgecolor="white",
        linewidth=0.3,
        label=name
    )
    left += vals

ax2.set_xlabel("Overlap genes by module")
ax2.set_title("Module burden", fontsize=13)

ax2.tick_params(axis="y", left=False, labelleft=False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Legends
handles = [
    plt.Line2D([0], [0], marker="o", color="w", label=k,
               markerfacecolor=v, markeredgecolor="black", markersize=8)
    for k, v in effect_colors.items()
]

leg1 = ax1.legend(
    handles=handles,
    title="Effect class",
    frameon=False,
    loc="lower right",
    fontsize=8,
    title_fontsize=9
)

ax2.legend(
    title="Gene module",
    frameon=False,
    bbox_to_anchor=(1.02, 1.0),
    loc="upper left",
    fontsize=7,
    title_fontsize=8
)

plt.tight_layout()

plt.savefig(os.path.join(FIG7, "7B_ranked_perturbation_landscape.tiff"), dpi=600, bbox_inches="tight")
plt.savefig(os.path.join(FIG7, "7B_ranked_perturbation_landscape.pdf"), bbox_inches="tight")

plt.close()

print("DONE 7B PLOT")
print(os.path.join(FIG7, "7B_ranked_perturbation_landscape.tiff"))