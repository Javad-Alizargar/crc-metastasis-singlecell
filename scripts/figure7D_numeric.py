import os
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

edges_file = os.path.join(FIG7, "7C_chord_edges.csv")
cand_file = os.path.join(FIG7, "7B_top50_candidates_cleaned.csv")

edges = pd.read_csv(edges_file)
cand = pd.read_csv(cand_file)

# ------------------------------------------------------------
# 1. Top drugs (same as 7C)
# ------------------------------------------------------------
top_drugs = cand.head(12)["clean_term"].tolist()

edges = edges[edges["drug_term"].isin(top_drugs)].copy()

# ------------------------------------------------------------
# 2. Gene recurrence
# ------------------------------------------------------------
gene_counts = edges["gene"].value_counts().reset_index()
gene_counts.columns = ["gene", "count"]

top_genes = gene_counts.head(30)["gene"].tolist()

edges = edges[edges["gene"].isin(top_genes)].copy()

# ------------------------------------------------------------
# 3. Node table
# ------------------------------------------------------------
nodes = []

# drug nodes
for d in top_drugs:
    score = cand[cand["clean_term"] == d]["final_rank_score"].values[0]
    effect = cand[cand["clean_term"] == d]["effect_class"].values[0]

    nodes.append({
        "id": d,
        "type": "drug",
        "size": score,
        "group": effect
    })

# gene nodes
for _, row in gene_counts.head(30).iterrows():
    nodes.append({
        "id": row["gene"],
        "type": "gene",
        "size": row["count"],
        "group": edges[edges["gene"] == row["gene"]]["module"].iloc[0]
    })

nodes = pd.DataFrame(nodes)

# ------------------------------------------------------------
# 4. Edge table
# ------------------------------------------------------------
edge_list = edges.copy()

edge_list = edge_list.rename(columns={
    "drug_term": "source",
    "gene": "target"
})

edge_list["weight"] = 1

edge_list = edge_list[["source", "target", "weight", "module"]]

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
nodes.to_csv(os.path.join(FIG7, "7D_network_nodes.csv"), index=False)
edge_list.to_csv(os.path.join(FIG7, "7D_network_edges.csv"), index=False)

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------
report = os.path.join(FIG7, "7D_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:

    f.write("FIGURE 7D NETWORK REPORT\n\n")

    f.write("Top drugs:\n")
    for d in top_drugs:
        f.write(d + "\n")

    f.write("\nTop genes:\n")
    for g in top_genes:
        f.write(g + "\n")

    f.write("\nGene hub counts:\n")
    f.write(gene_counts.head(20).to_string())

print("DONE 7D NUMERIC")
print(report)
print(open(report).read())