import os
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

cand_file = os.path.join(FIG7, "7B_top50_candidates_cleaned.csv")
overlap_file = os.path.join(FIG7, "7C_candidate_drug_overlap_genes_long.csv")
sig_file = os.path.join(
    OUT,
    "FIGURE_7_DRUG_REVERSAL_REBUILD",
    "7A_metastasis_signature_all_genes_annotated.csv"
)

cand = pd.read_csv(cand_file)
overlap = pd.read_csv(overlap_file)
sig = pd.read_csv(sig_file)

# ------------------------------------------------------------
# 1. Select top perturbations
# ------------------------------------------------------------
top_drugs = cand.head(12)["clean_term"].tolist()

sub = overlap[overlap["drug_term"].isin(top_drugs)].copy()

# ------------------------------------------------------------
# 2. Gene recurrence (importance)
# ------------------------------------------------------------
gene_counts = sub["gene"].value_counts().reset_index()
gene_counts.columns = ["gene", "count"]

top_genes = gene_counts.head(30)["gene"].tolist()

sub = sub[sub["gene"].isin(top_genes)].copy()

# ------------------------------------------------------------
# 3. Map modules
# ------------------------------------------------------------
gene_to_module = dict(zip(sig["gene"], sig["module"]))

sub["module"] = sub["gene"].map(gene_to_module).fillna("Other")

# ------------------------------------------------------------
# 4. Edge weight
# ------------------------------------------------------------
sub["weight"] = 1

# ------------------------------------------------------------
# 5. Save chord matrix (drug ↔ gene)
# ------------------------------------------------------------
matrix = pd.crosstab(sub["drug_term"], sub["gene"])

matrix.to_csv(os.path.join(FIG7, "7C_chord_matrix.csv"))

# ------------------------------------------------------------
# 6. Edge list (for advanced plotting)
# ------------------------------------------------------------
edges = sub[["drug_term", "gene", "module", "weight"]].copy()

edges.to_csv(os.path.join(FIG7, "7C_chord_edges.csv"), index=False)

# ------------------------------------------------------------
# 7. Node metadata
# ------------------------------------------------------------
nodes = []

# drug nodes
for d in top_drugs:
    nodes.append({
        "node": d,
        "type": "drug"
    })

# gene nodes
for g in top_genes:
    nodes.append({
        "node": g,
        "type": "gene",
        "module": gene_to_module.get(g, "Other")
    })

nodes = pd.DataFrame(nodes)

nodes.to_csv(os.path.join(FIG7, "7C_nodes.csv"), index=False)

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------
report = os.path.join(FIG7, "7C_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:

    f.write("FIGURE 7C CHORD NUMERIC REPORT\n\n")

    f.write("Top drugs:\n")
    for d in top_drugs:
        f.write(d + "\n")

    f.write("\nTop genes:\n")
    for g in top_genes:
        f.write(g + "\n")

    f.write("\nModule counts:\n")
    f.write(sub["module"].value_counts().to_string())

print("DONE 7C NUMERIC")
print(report)
print(open(report).read())