import os
import pandas as pd
import plotly.graph_objects as go

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

flow_file = os.path.join(FIG7, "7E_flow_drug_module_mechanism.csv")
df = pd.read_csv(flow_file)

# ----------------------------
# Clean drug names
# ----------------------------
def short_drug(x):
    x = str(x)
    x = x.replace("LJP006 ", "")
    x = x.replace("LJP005 ", "")
    x = x.replace("24H-", "")
    x = x.replace("3H-", "")
    x = x.replace("withaferin-a", "withaferin")
    x = x.replace("tanespimycin", "tanesp.")
    x = x.replace("geldanamycin", "geldan.")
    x = x.replace("puromycin", "puro.")
    x = x.replace("menadione", "menad.")
    x = x.replace("sanguinarine", "sang.")
    return x

df["drug_short"] = df["drug"].apply(short_drug)

# ----------------------------
# Nodes
# ----------------------------
drugs = list(df["drug_short"].unique())
modules = list(df["module"].unique())
mechanisms = list(df["mechanism"].unique())

labels = drugs + modules + mechanisms
label_to_index = {l: i for i, l in enumerate(labels)}

# ----------------------------
# Module colors (ONLY for flows)
# ----------------------------
module_colors = {
    "Stress_heatshock": "rgba(228,26,28,0.75)",        # red
    "Signaling_regulatory": "rgba(152,78,163,0.75)",   # purple
    "Other": "rgba(160,160,160,0.6)",                  # grey
    "Metabolic_lipid": "rgba(77,175,74,0.75)",         # green
    "Antigen_presentation_myeloid": "rgba(166,86,40,0.75)"
}

# ----------------------------
# Build links
# ----------------------------
sources, targets, values, colors = [], [], [], []

# Drug → Module
dm = df.groupby(["drug_short", "module"])["weight"].sum().reset_index()

for _, row in dm.iterrows():
    sources.append(label_to_index[row["drug_short"]])
    targets.append(label_to_index[row["module"]])
    values.append(row["weight"])
    colors.append(module_colors.get(row["module"], "rgba(150,150,150,0.5)"))

# Module → Mechanism
mm = df.groupby(["module", "mechanism"])["weight"].sum().reset_index()

for _, row in mm.iterrows():
    sources.append(label_to_index[row["module"]])
    targets.append(label_to_index[row["mechanism"]])
    values.append(row["weight"])
    colors.append(module_colors.get(row["module"], "rgba(150,150,150,0.5)"))

# ----------------------------
# Sankey plot (clean style)
# ----------------------------
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=34,
        thickness=28,
        line=dict(color="black", width=0.8),
        label=labels,
        color=[
            "rgba(245,245,245,0.95)" if x in drugs else
            "rgba(255,255,255,0.98)" if x in modules else
            "rgba(245,245,245,0.95)"
            for x in labels
        ]
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=colors
    )
)])

# ----------------------------
# Layout (BIG readable text)
# ----------------------------
fig.update_layout(
    width=2100,
    height=1350,
    font=dict(size=24, color="black", family="Arial"),
    margin=dict(l=40, r=40, t=40, b=40),
    paper_bgcolor="white",
    plot_bgcolor="white"
)

# ----------------------------
# Save (HIGH QUALITY)
# ----------------------------
fig.write_image(
    os.path.join(FIG7, "7E_mechanism_sankey_FINAL.png"),
    scale=3
)

fig.write_image(
    os.path.join(FIG7, "7E_mechanism_sankey_FINAL.pdf")
)

print("DONE FINAL CLEAN 7E")
print(os.path.join(FIG7, "7E_mechanism_sankey_FINAL.png"))