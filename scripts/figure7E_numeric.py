import os
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

edges = pd.read_csv(os.path.join(FIG7, "7C_chord_edges.csv"))

# ------------------------------------------------------------
# Top drugs + genes (same filtering as before)
# ------------------------------------------------------------
cand = pd.read_csv(os.path.join(FIG7, "7B_top50_candidates_cleaned.csv"))
top_drugs = cand.head(12)["clean_term"].tolist()

edges = edges[edges["drug_term"].isin(top_drugs)].copy()

gene_counts = edges["gene"].value_counts().reset_index()
gene_counts.columns = ["gene", "count"]

top_genes = gene_counts.head(30)["gene"].tolist()
edges = edges[edges["gene"].isin(top_genes)].copy()

# ------------------------------------------------------------
# Module → mechanism mapping
# ------------------------------------------------------------
module_to_mechanism = {
    "Stress_heatshock": "Protein folding / stress response",
    "Signaling_regulatory": "Signal transduction / transcriptional control",
    "Cytotoxic_NK_T": "Cytotoxic immune activation",
    "Cytokine_inflammatory": "Cytokine signaling / inflammation",
    "Metabolic_lipid": "Metabolic reprogramming",
    "Antigen_presentation_myeloid": "Myeloid activation / antigen presentation",
    "Other": "Other / mixed"
}

edges["mechanism"] = edges["module"].map(module_to_mechanism).fillna("Other / mixed")

# ------------------------------------------------------------
# Build flows
# ------------------------------------------------------------
flow_rows = []

for _, row in edges.iterrows():
    flow_rows.append({
        "drug": row["drug_term"],
        "module": row["module"],
        "mechanism": row["mechanism"],
        "weight": 1
    })

flow = pd.DataFrame(flow_rows)

# Aggregate
flow_agg = (
    flow.groupby(["drug", "module", "mechanism"])["weight"]
    .sum()
    .reset_index()
)

# Also module→mechanism summary
module_mech = (
    flow.groupby(["module", "mechanism"])["weight"]
    .sum()
    .reset_index()
)

# Save
flow_agg.to_csv(os.path.join(FIG7, "7E_flow_drug_module_mechanism.csv"), index=False)
module_mech.to_csv(os.path.join(FIG7, "7E_module_mechanism_summary.csv"), index=False)

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------
report = os.path.join(FIG7, "7E_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:

    f.write("FIGURE 7E MECHANISM FLOW REPORT\n\n")

    f.write("Module → mechanism counts:\n")
    f.write(module_mech.to_string())

    f.write("\n\nTop contributing flows:\n")
    f.write(flow_agg.sort_values("weight", ascending=False).head(20).to_string())

print("DONE 7E NUMERIC")
print(report)