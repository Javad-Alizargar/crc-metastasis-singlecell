import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL")
os.makedirs(FIG7, exist_ok=True)

# ------------------------------------------------------------
# Load DE results
# ------------------------------------------------------------
de_file = os.path.join(
    OUT,
    "FIGURE_5_METASTASIS_MECHANISM",
    "Figure5A_primary_vs_metastasis_DE_all_genes.csv"
)

de = pd.read_csv(de_file)

# ------------------------------------------------------------
# Metastasis signature
# ------------------------------------------------------------
meta_genes = (
    de[de["direction"] == "metastasis_up"]
    .sort_values("log2FC", ascending=False)
    .head(200)
)

meta_genes["weight"] = meta_genes["log2FC"]

meta_genes.to_csv(os.path.join(FIG7, "7A_metastasis_signature_genes.csv"), index=False)

# ------------------------------------------------------------
# MOCK DRUG DATABASE (replace later if needed)
# ------------------------------------------------------------
# Each drug = genes it DOWNREGULATES
drug_db = {
    "Drug_A_JAK_inhibitor": ["IL12RB2", "CCL5", "XCL1", "XCL2"],
    "Drug_B_HSP_inhibitor": ["HSPA1A", "HSPA6", "HSPH1", "HSPD1"],
    "Drug_C_NFkB_inhibitor": ["TNFSF13", "TNFSF14", "CCL5"],
    "Drug_D_KLR_pathway_block": ["KLRF1", "KLRC1", "KLRC4"],
    "Drug_E_metabolic": ["CYP51A1", "SLC4A4"]
}

# ------------------------------------------------------------
# Drug reversal scoring
# ------------------------------------------------------------
results = []

for drug, targets in drug_db.items():

    overlap = meta_genes[meta_genes["gene"].isin(targets)]

    if len(overlap) == 0:
        continue

    score = overlap["weight"].sum()
    n_overlap = len(overlap)

    results.append({
        "drug": drug,
        "overlap_genes": n_overlap,
        "reversal_score": score,
        "genes": ",".join(overlap["gene"].tolist())
    })

drug_df = pd.DataFrame(results)

drug_df = drug_df.sort_values("reversal_score", ascending=False)

drug_df.to_csv(os.path.join(FIG7, "7B_drug_ranking.csv"), index=False)

# ------------------------------------------------------------
# 7B — Barplot
# ------------------------------------------------------------
plt.figure(figsize=(6,4))

plt.barh(
    drug_df["drug"],
    drug_df["reversal_score"]
)

plt.xlabel("Reversal score")
plt.ylabel("")
plt.title("")
plt.gca().invert_yaxis()

plt.tight_layout()

plt.savefig(os.path.join(FIG7, "7B_drug_ranking.tiff"), dpi=600)
plt.savefig(os.path.join(FIG7, "7B_drug_ranking.pdf"))

plt.close()

# ------------------------------------------------------------
# 7C — Heatmap (drug vs gene)
# ------------------------------------------------------------
genes = meta_genes["gene"].head(30).tolist()

heat = pd.DataFrame(0, index=drug_df["drug"], columns=genes)

for drug, targets in drug_db.items():
    for g in targets:
        if g in heat.columns and drug in heat.index:
            heat.loc[drug, g] = 1

plt.figure(figsize=(10,4))
plt.imshow(heat, aspect="auto")

plt.xticks(range(len(genes)), genes, rotation=90)
plt.yticks(range(len(heat.index)), heat.index)

plt.colorbar(label="Target")

plt.tight_layout()

plt.savefig(os.path.join(FIG7, "7C_drug_gene_heatmap.tiff"), dpi=600)
plt.savefig(os.path.join(FIG7, "7C_drug_gene_heatmap.pdf"))

plt.close()

# ------------------------------------------------------------
# 7D — Mechanism grouping
# ------------------------------------------------------------
drug_df["category"] = drug_df["drug"].apply(lambda x:
    "Stress" if "HSP" in x else
    "Immune" if "JAK" in x or "NFkB" in x else
    "Cytotoxic" if "KLR" in x else
    "Metabolic"
)

cat_summary = drug_df.groupby("category")["reversal_score"].mean()

cat_summary.to_csv(os.path.join(FIG7, "7D_category_summary.csv"))

# ------------------------------------------------------------
# 7E — Top target genes
# ------------------------------------------------------------
top_targets = []

for _, row in drug_df.iterrows():
    genes = row["genes"].split(",")
    for g in genes:
        top_targets.append(g)

top_targets = pd.Series(top_targets).value_counts().head(20)

top_targets.to_csv(os.path.join(FIG7, "7E_top_targets.csv"))

print("DONE FIGURE 7")
print(drug_df)