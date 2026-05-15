import os
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG9 = os.path.join(OUT, "FIGURE_9_CONCEPTUAL_MODEL")
os.makedirs(FIG9, exist_ok=True)

# ------------------------------------------------------------
# Evidence matrix: manually curated from your completed figures
# Score:
# 0 = no direct evidence
# 1 = weak / indirect evidence
# 2 = moderate evidence
# 3 = strong evidence
# ------------------------------------------------------------

rows = [
    {
        "claim": "Stress-response axis",
        "short_label": "Stress axis",
        "DE": 3,
        "Pathway_enrichment": 3,
        "UMAP_score": 3,
        "Pseudotime": 1,
        "Perturbation_network": 3,
        "Key_evidence": "HSPA1A, HSPA6, HSPH1, HSPD1, DNAJB1, HSP90AA1; enriched protein folding and unfolded protein response; dominant perturbation/network module."
    },
    {
        "claim": "Immune activation in metastasis",
        "short_label": "Immune activation",
        "DE": 3,
        "Pathway_enrichment": 3,
        "UMAP_score": 3,
        "Pseudotime": 2,
        "Perturbation_network": 1,
        "Key_evidence": "XCL1, XCL2, CCL5, KLRF1, KLRC1, IL12RB2; enriched cytokine production, NK cytotoxicity, inflammatory response."
    },
    {
        "claim": "Myeloid/macrophage late immune state",
        "short_label": "Myeloid late state",
        "DE": 2,
        "Pathway_enrichment": 2,
        "UMAP_score": 1,
        "Pseudotime": 3,
        "Perturbation_network": 1,
        "Key_evidence": "Late pseudotime Q4 enriched for macrophage and myeloid cells; myeloid mean pseudotime highest, followed by macrophage."
    },
    {
        "claim": "Plasma-cell early immune state",
        "short_label": "Plasma early state",
        "DE": 1,
        "Pathway_enrichment": 1,
        "UMAP_score": 1,
        "Pseudotime": 3,
        "Perturbation_network": 0,
        "Key_evidence": "Early pseudotime Q1 enriched for plasma cells; plasma cells show strongest primary-to-metastasis pseudotime shift."
    },
    {
        "claim": "Primary epithelial/structural program",
        "short_label": "Primary structure",
        "DE": 3,
        "Pathway_enrichment": 3,
        "UMAP_score": 3,
        "Pseudotime": 1,
        "Perturbation_network": 0,
        "Key_evidence": "IGFBP6, CAV1, CAV2, COL6A3, KRT8, TM4SF1; enriched focal adhesion, ECM-receptor interaction, EMT and structural programs."
    },
    {
        "claim": "Limited tumor-epithelial trajectory",
        "short_label": "Tumor stability",
        "DE": 2,
        "Pathway_enrichment": 1,
        "UMAP_score": 1,
        "Pseudotime": 3,
        "Perturbation_network": 0,
        "Key_evidence": "Tumor epithelial primary and metastasis pseudotime means are close; primary and metastatic tumor cells overlap substantially."
    },
    {
        "claim": "Perturbation hits mostly mimic stress program",
        "short_label": "Drug mimicry",
        "DE": 1,
        "Pathway_enrichment": 1,
        "UMAP_score": 0,
        "Pseudotime": 0,
        "Perturbation_network": 3,
        "Key_evidence": "Most top perturbation hits are UP signatures classified as potential risk/mimic; stress heat-shock genes dominate overlap network."
    },
    {
        "claim": "Candidate therapeutic target axis",
        "short_label": "Target axis",
        "DE": 2,
        "Pathway_enrichment": 2,
        "UMAP_score": 2,
        "Pseudotime": 1,
        "Perturbation_network": 3,
        "Key_evidence": "HSPA1A, HSPA6, HSPH1, DNAJB1, BAG3, HMOX1 repeatedly appear as hubs across DE, pathways, perturbation chord, and network analysis."
    }
]

evidence = pd.DataFrame(rows)

# Long format for ggplot heatmap
long = evidence.melt(
    id_vars=["claim", "short_label", "Key_evidence"],
    value_vars=["DE", "Pathway_enrichment", "UMAP_score", "Pseudotime", "Perturbation_network"],
    var_name="Evidence_layer",
    value_name="Strength"
)

# Friendly labels
layer_labels = {
    "DE": "Differential expression",
    "Pathway_enrichment": "Pathway enrichment",
    "UMAP_score": "Score overlays",
    "Pseudotime": "Pseudotime",
    "Perturbation_network": "Perturbation network"
}
long["Evidence_layer_label"] = long["Evidence_layer"].map(layer_labels)

strength_labels = {
    0: "None",
    1: "Weak",
    2: "Moderate",
    3: "Strong"
}
long["Strength_label"] = long["Strength"].map(strength_labels)

# Save
evidence.to_csv(os.path.join(FIG9, "Figure9B_evidence_matrix_wide.csv"), index=False)
long.to_csv(os.path.join(FIG9, "Figure9B_evidence_matrix_long.csv"), index=False)

# TXT report
report = os.path.join(FIG9, "Figure9B_EVIDENCE_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:
    f.write("FIGURE 9B EVIDENCE MATRIX NUMERIC REPORT\n\n")

    f.write("1. Evidence strength coding\n")
    f.write("0 = none\n")
    f.write("1 = weak / indirect\n")
    f.write("2 = moderate\n")
    f.write("3 = strong\n\n")

    f.write("2. Wide evidence matrix\n")
    f.write(evidence.to_string(index=False))
    f.write("\n\n")

    f.write("3. Long evidence matrix preview\n")
    f.write(long.head(30).to_string(index=False))
    f.write("\n\n")

    f.write("4. Files generated\n")
    for x in sorted(os.listdir(FIG9)):
        f.write(x + "\n")

print("DONE FIGURE 9B NUMERIC")
print(report)
print(open(report, encoding="utf-8").read())