import os
import textwrap
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

nodes = pd.read_csv(os.path.join(FIG7, "7D_network_nodes.csv"))
edges = pd.read_csv(os.path.join(FIG7, "7D_network_edges.csv"))

G = nx.Graph()

for _, row in nodes.iterrows():
    G.add_node(
        row["id"],
        node_type=row["type"],
        group=row["group"],
        size=row["size"]
    )

for _, row in edges.iterrows():
    G.add_edge(
        row["source"],
        row["target"],
        weight=row["weight"],
        module=row["module"]
    )

module_colors = {
    "Stress_heatshock": "#E41A1C",
    "Signaling_regulatory": "#984EA3",
    "Other": "#BDBDBD",
    "Metabolic_lipid": "#4DAF4A",
    "Antigen_presentation_myeloid": "#A65628",
    "Cytotoxic_NK_T": "#377EB8",
    "Cytokine_inflammatory": "#FF7F00",
    "potential_risk": "#D73027",
    "candidate_reversal": "#1A9850",
    "association": "#6A3D9A"
}

# More spread-out layout
pos = nx.spring_layout(
    G,
    k=1.15,
    iterations=500,
    seed=42,
    weight="weight"
)

drug_nodes = [n for n, d in G.nodes(data=True) if d["node_type"] == "drug"]
gene_nodes = [n for n, d in G.nodes(data=True) if d["node_type"] == "gene"]

drug_max = max(nodes[nodes["type"] == "drug"]["size"])
gene_max = max(nodes[nodes["type"] == "gene"]["size"])

drug_sizes = [
    700 + (G.nodes[n]["size"] / drug_max) * 1400
    for n in drug_nodes
]

gene_sizes = [
    650 + (G.nodes[n]["size"] / gene_max) * 1600
    for n in gene_nodes
]

drug_colors = [module_colors.get(G.nodes[n]["group"], "#999999") for n in drug_nodes]
gene_colors = [module_colors.get(G.nodes[n]["group"], "#999999") for n in gene_nodes]

edge_colors = [
    module_colors.get(d.get("module", "Other"), "#BDBDBD")
    for _, _, d in G.edges(data=True)
]

plt.figure(figsize=(18, 15))

# Edges
nx.draw_networkx_edges(
    G,
    pos,
    edge_color=edge_colors,
    alpha=0.22,
    width=1.0
)

# Drug nodes = squares
nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=drug_nodes,
    node_size=drug_sizes,
    node_color=drug_colors,
    edgecolors="black",
    linewidths=1.1,
    node_shape="s",
    alpha=0.92
)

# Gene nodes = circles
nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=gene_nodes,
    node_size=gene_sizes,
    node_color=gene_colors,
    edgecolors="black",
    linewidths=0.6,
    node_shape="o",
    alpha=0.95
)

# Drug labels: wrapped, outside/near squares
drug_labels = {
    n: "\n".join(textwrap.wrap(n, width=18))
    for n in drug_nodes
}

nx.draw_networkx_labels(
    G,
    pos,
    labels=drug_labels,
    font_size=6.5,
    font_color="black",
    verticalalignment="center"
)

# Gene labels: all gene circles labeled
gene_labels = {n: n for n in gene_nodes}

nx.draw_networkx_labels(
    G,
    pos,
    labels=gene_labels,
    font_size=7,
    font_color="black",
    font_weight="bold",
    verticalalignment="center"
)

plt.axis("off")

# Legend OUTSIDE right side, no overlap
legend_items = [
    ("Stress / heat shock genes", "#E41A1C"),
    ("Signaling regulatory genes", "#984EA3"),
    ("Other genes", "#BDBDBD"),
    ("Metabolic / lipid genes", "#4DAF4A"),
    ("Antigen / myeloid genes", "#A65628"),
    ("Potential risk / mimic perturbations", "#D73027"),
    ("Candidate reversal perturbations", "#1A9850"),
    ("Association", "#6A3D9A"),
]

for label, color in legend_items:
    plt.scatter([], [], c=color, s=120, edgecolors="black", label=label)

plt.legend(
    frameon=False,
    bbox_to_anchor=(1.02, 0.5),
    loc="center left",
    fontsize=8
)

plt.tight_layout(rect=[0, 0, 0.82, 1])

plt.savefig(
    os.path.join(FIG7, "7D_drug_gene_network_CLEAN.tiff"),
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(FIG7, "7D_drug_gene_network_CLEAN.pdf"),
    bbox_inches="tight"
)

plt.close()

print("DONE CLEAN 7D network")
print(os.path.join(FIG7, "7D_drug_gene_network_CLEAN.tiff"))