import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_ADVANCED")
os.makedirs(FIG7, exist_ok=True)

de_file = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
de = pd.read_csv(de_file)

meta_genes = (
    de[de["direction"] == "metastasis_up"]
    .sort_values("log2FC", ascending=False)["gene"]
    .head(300)
    .dropna()
    .astype(str)
    .tolist()
)

primary_genes = (
    de[de["direction"] == "primary_up"]
    .sort_values("log2FC", ascending=True)["gene"]
    .head(300)
    .dropna()
    .astype(str)
    .tolist()
)

pd.DataFrame({"metastasis_up_genes": meta_genes}).to_csv(
    os.path.join(FIG7, "7A_metastasis_signature_top300.csv"),
    index=False
)

pd.DataFrame({"primary_up_genes": primary_genes}).to_csv(
    os.path.join(FIG7, "7A_primary_signature_top300.csv"),
    index=False
)

# ------------------------------------------------------------
# Find available drug libraries automatically
# ------------------------------------------------------------
libs = gp.get_library_name(organism="human")

wanted_keywords = [
    "LINCS",
    "Chem",
    "Drug",
    "Perturb",
    "DSigDB",
    "CMap"
]

drug_libs = [
    x for x in libs
    if any(k.lower() in x.lower() for k in wanted_keywords)
]

pd.DataFrame({"available_drug_like_libraries": drug_libs}).to_csv(
    os.path.join(FIG7, "available_drug_like_libraries.csv"),
    index=False
)

print("Available drug-like libraries:")
for x in drug_libs:
    print(" -", x)

# Prefer true reversal libraries if present
preferred_libs = []
for target in [
    "LINCS_L1000_Chem_Pert_down",
    "LINCS_L1000_Chem_Pert_up",
    "DSigDB",
    "Drug_Perturbations_from_GEO_2014"
]:
    if target in libs:
        preferred_libs.append(target)

# Fallback: use any drug-like libraries
if len(preferred_libs) == 0:
    preferred_libs = drug_libs[:5]

print("\nUsing libraries:")
for x in preferred_libs:
    print(" -", x)

# ------------------------------------------------------------
# Enrichment helper
# ------------------------------------------------------------
def run_enrichr(gene_list, label):
    all_results = []

    for lib in preferred_libs:
        print("Running", label, lib)

        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism="human",
                outdir=None,
                cutoff=1.0
            )

            res = enr.results.copy()
            if res is None or res.shape[0] == 0:
                continue

            res["query_signature"] = label
            res["library"] = lib

            all_results.append(res)

        except Exception as e:
            print("FAILED", lib, e)

    if len(all_results) == 0:
        return pd.DataFrame()

    out = pd.concat(all_results, ignore_index=True)
    out["Adjusted P-value"] = pd.to_numeric(out["Adjusted P-value"], errors="coerce")
    out["Combined Score"] = pd.to_numeric(out["Combined Score"], errors="coerce")
    out = out.sort_values(["Adjusted P-value", "Combined Score"], ascending=[True, False])

    return out

meta_drugs = run_enrichr(meta_genes, "metastasis_up")
primary_drugs = run_enrichr(primary_genes, "primary_up")

meta_drugs.to_csv(os.path.join(FIG7, "7B_metastasis_up_drug_enrichment_all.csv"), index=False)
primary_drugs.to_csv(os.path.join(FIG7, "7B_primary_up_drug_enrichment_all.csv"), index=False)

# ------------------------------------------------------------
# Define reversal candidate logic
# Disease metastasis-up genes enriched in drug DOWN signatures
# are stronger reversal candidates.
# ------------------------------------------------------------
def classify_reversal(row):
    lib = str(row.get("library", "")).lower()
    term = str(row.get("Term", "")).lower()

    if "down" in lib or "down" in term:
        return "strong_reversal_candidate"
    if "up" in lib or "up" in term:
        return "opposite_direction_or_risk"
    return "drug_association"

if meta_drugs.shape[0] > 0:
    meta_drugs["reversal_class"] = meta_drugs.apply(classify_reversal, axis=1)
    meta_drugs["minus_log10_adj_p"] = -np.log10(meta_drugs["Adjusted P-value"].clip(lower=1e-300))
    meta_drugs["rank_score"] = meta_drugs["minus_log10_adj_p"] * np.log1p(meta_drugs["Combined Score"])

    candidates = meta_drugs.sort_values("rank_score", ascending=False).head(50).copy()
else:
    candidates = pd.DataFrame()

candidates.to_csv(os.path.join(FIG7, "7B_top50_candidate_drug_reversal.csv"), index=False)

# ------------------------------------------------------------
# Extract overlap genes
# ------------------------------------------------------------
def split_genes(x):
    if pd.isna(x):
        return []
    x = str(x)
    parts = re.split(r"[;,/]", x)
    parts = [p.strip() for p in parts if p.strip()]
    return parts

if candidates.shape[0] > 0:
    candidate_gene_rows = []
    for _, row in candidates.iterrows():
        genes = split_genes(row.get("Genes", ""))
        for g in genes:
            candidate_gene_rows.append({
                "drug_term": row["Term"],
                "library": row["library"],
                "gene": g,
                "adjusted_p": row["Adjusted P-value"],
                "combined_score": row["Combined Score"],
                "rank_score": row["rank_score"],
                "reversal_class": row["reversal_class"]
            })

    candidate_gene_df = pd.DataFrame(candidate_gene_rows)
else:
    candidate_gene_df = pd.DataFrame()

candidate_gene_df.to_csv(os.path.join(FIG7, "7C_candidate_drug_overlap_genes_long.csv"), index=False)

# ------------------------------------------------------------
# PANEL 7B: ranked drug candidates
# ------------------------------------------------------------
if candidates.shape[0] > 0:
    top_plot = candidates.head(25).copy()
    top_plot = top_plot.sort_values("rank_score", ascending=True)

    plt.figure(figsize=(10, 8))
    plt.barh(top_plot["Term"], top_plot["rank_score"], color="#6A3D9A")
    plt.xlabel("Drug reversal rank score")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG7, "7B_ranked_drug_candidates.tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG7, "7B_ranked_drug_candidates.pdf"), bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# PANEL 7C: significance vs combined score bubble plot
# ------------------------------------------------------------
if candidates.shape[0] > 0:
    plot = candidates.head(40).copy()
    plot["neglogp"] = -np.log10(plot["Adjusted P-value"].clip(lower=1e-300))

    plt.figure(figsize=(8, 6))
    plt.scatter(
        plot["Combined Score"],
        plot["neglogp"],
        s=np.clip(plot["Overlap"].astype(str).str.extract(r"(\d+)")[0].fillna(5).astype(float) * 20, 50, 500),
        alpha=0.65,
        c="#D73027"
    )

    for _, row in plot.head(12).iterrows():
        plt.text(row["Combined Score"], row["neglogp"], str(row["Term"])[:35], fontsize=7)

    plt.xlabel("Combined score")
    plt.ylabel("-log10 adjusted p-value")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG7, "7C_drug_significance_bubble.tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG7, "7C_drug_significance_bubble.pdf"), bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# PANEL 7D: drug-term vs overlap-gene heatmap
# ------------------------------------------------------------
if candidate_gene_df.shape[0] > 0:
    top_terms = candidates.head(20)["Term"].tolist()
    top_genes = (
        candidate_gene_df[candidate_gene_df["drug_term"].isin(top_terms)]["gene"]
        .value_counts()
        .head(35)
        .index
        .tolist()
    )

    heat = pd.DataFrame(0, index=top_terms, columns=top_genes)

    for _, row in candidate_gene_df.iterrows():
        if row["drug_term"] in heat.index and row["gene"] in heat.columns:
            heat.loc[row["drug_term"], row["gene"]] = 1

    heat.to_csv(os.path.join(FIG7, "7D_drug_gene_overlap_matrix.csv"))

    plt.figure(figsize=(13, 8))
    plt.imshow(heat, aspect="auto", cmap="Reds")
    plt.xticks(range(len(heat.columns)), heat.columns, rotation=90, fontsize=7)
    plt.yticks(range(len(heat.index)), [x[:45] for x in heat.index], fontsize=7)
    plt.colorbar(label="Overlap")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG7, "7D_drug_gene_overlap_heatmap.tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG7, "7D_drug_gene_overlap_heatmap.pdf"), bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# PANEL 7E: recurrent target genes
# ------------------------------------------------------------
if candidate_gene_df.shape[0] > 0:
    target_counts = candidate_gene_df["gene"].value_counts().head(30).reset_index()
    target_counts.columns = ["gene", "drug_candidate_count"]
    target_counts.to_csv(os.path.join(FIG7, "7E_recurrent_drug_target_genes.csv"), index=False)

    plot = target_counts.sort_values("drug_candidate_count", ascending=True)

    plt.figure(figsize=(7, 7))
    plt.barh(plot["gene"], plot["drug_candidate_count"], color="#1F78B4")
    plt.xlabel("Number of candidate drug signatures")
    plt.ylabel("")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG7, "7E_recurrent_target_genes.tiff"), dpi=600, bbox_inches="tight")
    plt.savefig(os.path.join(FIG7, "7E_recurrent_target_genes.pdf"), bbox_inches="tight")
    plt.close()

# ------------------------------------------------------------
# REPORT
# ------------------------------------------------------------
report_path = os.path.join(FIG7, "FIGURE7_DRUG_REVERSAL_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 7 ADVANCED DRUG REVERSAL REPORT\n\n")

    f.write("1. Libraries used\n")
    f.write("\n".join(preferred_libs))
    f.write("\n\n")

    f.write("2. Metastasis signature genes used\n")
    f.write(", ".join(meta_genes[:100]))
    f.write("\n\n")

    f.write("3. Top candidate drug/reversal terms\n")
    if candidates.shape[0] > 0:
        cols = ["Term", "library", "Adjusted P-value", "Combined Score", "Overlap", "Genes", "reversal_class", "rank_score"]
        cols = [c for c in cols if c in candidates.columns]
        f.write(candidates[cols].head(30).to_string(index=False))
    else:
        f.write("No candidates found.")
    f.write("\n\n")

    f.write("4. Recurrent overlap genes\n")
    if candidate_gene_df.shape[0] > 0:
        f.write(candidate_gene_df["gene"].value_counts().head(30).to_string())
    else:
        f.write("No overlap genes.")
    f.write("\n\n")

    f.write("5. Files generated\n")
    for x in sorted(os.listdir(FIG7)):
        f.write(x + "\n")

print("DONE ADVANCED FIGURE 7")
print(report_path)
print(open(report_path, encoding="utf-8").read())