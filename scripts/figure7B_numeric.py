import os
import re
import pandas as pd
import numpy as np

OUT = r"D:\CRC_META_FULL_SCVI"
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

cand_file = os.path.join(FIG7, "7B_top50_candidate_drug_reversal.csv")
gene_overlap_file = os.path.join(FIG7, "7C_candidate_drug_overlap_genes_long.csv")
sig_file = os.path.join(
    OUT,
    "FIGURE_7_DRUG_REVERSAL_REBUILD",
    "7A_metastasis_signature_all_genes_annotated.csv"
)

cand = pd.read_csv(cand_file)
overlap = pd.read_csv(gene_overlap_file)
sig = pd.read_csv(sig_file)

# ------------------------------------------------------------
# 1. Clean names
# ------------------------------------------------------------
def clean_term(x):
    x = str(x)
    x = re.sub(r"\(.*?\)", "", x)
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x).strip()
    return x

cand["clean_term"] = cand["Term"].apply(clean_term)

# ------------------------------------------------------------
# 2. Direction classification (CRITICAL)
# ------------------------------------------------------------
def classify_direction(row):
    txt = (str(row["Term"]) + " " + str(row["library"])).lower()
    if "down" in txt:
        return "DOWN"
    elif "up" in txt:
        return "UP"
    else:
        return "UNKNOWN"

cand["direction"] = cand.apply(classify_direction, axis=1)

# ------------------------------------------------------------
# 3. Reversal logic
# metastasis-up genes → DOWN perturbation = reversal
# ------------------------------------------------------------
def classify_effect(row):
    if row["direction"] == "DOWN":
        return "candidate_reversal"
    elif row["direction"] == "UP":
        return "potential_risk"
    else:
        return "association"

cand["effect_class"] = cand.apply(classify_effect, axis=1)

# ------------------------------------------------------------
# 4. Overlap gene parsing
# ------------------------------------------------------------
def split_genes(x):
    if pd.isna(x):
        return []
    parts = re.split(r"[;,/]", str(x))
    return [p.strip() for p in parts if p.strip()]

cand["gene_list"] = cand["Genes"].apply(split_genes)

# ------------------------------------------------------------
# 5. Map gene → module (from 7A)
# ------------------------------------------------------------
gene_to_module = dict(zip(sig["gene"], sig["module"]))

def count_modules(glist):
    counts = {}
    for g in glist:
        m = gene_to_module.get(g, "Other")
        counts[m] = counts.get(m, 0) + 1
    return counts

module_counts = cand["gene_list"].apply(count_modules)

# expand module counts into columns
all_modules = sorted(set(sig["module"].unique()))
for m in all_modules:
    cand[f"module_{m}"] = module_counts.apply(lambda x: x.get(m, 0))

# ------------------------------------------------------------
# 6. Derived scores
# ------------------------------------------------------------
cand["n_overlap"] = cand["gene_list"].apply(len)

cand["minus_log10_padj"] = -np.log10(
    cand["Adjusted P-value"].replace(0, 1e-300)
)

cand["final_rank_score"] = (
    cand["minus_log10_padj"] *
    np.log1p(cand["Combined Score"]) *
    np.log1p(cand["n_overlap"])
)

# ------------------------------------------------------------
# 7. Top selections
# ------------------------------------------------------------
cand_sorted = cand.sort_values("final_rank_score", ascending=False)

top20 = cand_sorted.head(20).copy()
top50 = cand_sorted.head(50).copy()

# ------------------------------------------------------------
# 8. Summary stats
# ------------------------------------------------------------
effect_counts = cand["effect_class"].value_counts()

direction_counts = cand["direction"].value_counts()

module_totals = {}
for m in all_modules:
    module_totals[m] = cand[f"module_{m}"].sum()

module_totals = pd.Series(module_totals).sort_values(ascending=False)

# ------------------------------------------------------------
# 9. Save CSVs
# ------------------------------------------------------------
cand_sorted.to_csv(os.path.join(FIG7, "7B_candidates_ranked_full.csv"), index=False)
top20.to_csv(os.path.join(FIG7, "7B_top20_candidates.csv"), index=False)
top50.to_csv(os.path.join(FIG7, "7B_top50_candidates_cleaned.csv"), index=False)

# ------------------------------------------------------------
# 10. TXT REPORT
# ------------------------------------------------------------
report = os.path.join(FIG7, "7B_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:

    f.write("FIGURE 7B DRUG PERTURBATION NUMERIC REPORT\n\n")

    f.write("1. Input files\n")
    f.write(cand_file + "\n\n")

    f.write("2. Effect classification\n")
    f.write(effect_counts.to_string())
    f.write("\n\n")

    f.write("3. Direction counts\n")
    f.write(direction_counts.to_string())
    f.write("\n\n")

    f.write("4. Top 20 candidates\n")
    cols = [
        "clean_term",
        "library",
        "direction",
        "effect_class",
        "Adjusted P-value",
        "Combined Score",
        "n_overlap",
        "final_rank_score"
    ]
    f.write(top20[cols].to_string(index=False))
    f.write("\n\n")

    f.write("5. Module contributions across all candidates\n")
    f.write(module_totals.to_string())
    f.write("\n\n")

    f.write("6. Example top candidate gene overlaps\n")
    for i, row in top20.head(5).iterrows():
        f.write("\n---\n")
        f.write(row["clean_term"] + "\n")
        f.write("Genes: " + ", ".join(row["gene_list"][:20]) + "\n")

print("DONE 7B NUMERIC")
print(report)
print(open(report, encoding="utf-8").read())