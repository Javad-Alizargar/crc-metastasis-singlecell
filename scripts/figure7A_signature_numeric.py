import os
import re
import numpy as np
import pandas as pd

OUT = r"D:\CRC_META_FULL_SCVI"
FIG5 = os.path.join(OUT, "FIGURE_5_METASTASIS_MECHANISM")
FIG7 = os.path.join(OUT, "FIGURE_7_DRUG_REVERSAL_REBUILD")
os.makedirs(FIG7, exist_ok=True)

de_file = os.path.join(FIG5, "Figure5A_primary_vs_metastasis_DE_all_genes.csv")
pathway_file = os.path.join(FIG5, "Figure5D_metastasis_up_top20.csv")

de = pd.read_csv(de_file)

meta = de[de["direction"] == "metastasis_up"].copy()
meta = meta.sort_values(["log2FC", "score"], ascending=[False, False]).reset_index(drop=True)

# ------------------------------------------------------------
# Gene module dictionaries based on your actual Figure 5 biology
# ------------------------------------------------------------
modules = {
    "Stress_heatshock": [
        "HSPA1A", "HSPA1B", "HSPA2", "HSPA5", "HSPA6", "HSPA8",
        "HSPH1", "HSPD1", "HSPE1", "HSP90AA1", "DNAJB1", "DNAJA4",
        "DNAJB4", "BAG3", "FKBP4", "HMOX1"
    ],
    "Cytotoxic_NK_T": [
        "KLRF1", "KLRC1", "KLRC2", "KLRC4", "KLRD1", "KLRB1",
        "GNLY", "GZMA", "GZMK", "PRF1", "CD160", "CD69",
        "TRGC1", "TRGV10", "TRBC1", "CD8A"
    ],
    "Cytokine_inflammatory": [
        "XCL1", "XCL2", "CCL5", "CCL4", "TNF", "TNFSF13",
        "TNFSF14", "IFNG", "IL12RB2", "IFNGR1", "CSF1", "OSM"
    ],
    "Antigen_presentation_myeloid": [
        "HLA-DRB5", "IFI30", "FCN1", "ITGAM", "FCER1G", "CEBPD",
        "CD300A", "C5AR1", "CD163", "MS4A6A", "TYROBP"
    ],
    "Metabolic_lipid": [
        "CYP51A1", "SLC4A4", "GPAT3", "AGPAT4", "BCO2", "ABCB1",
        "CTH", "INSIG1", "ACSL1", "HIBCH", "PLIN2"
    ],
    "Signaling_regulatory": [
        "YES1", "VAV3", "PLCG2", "PLCB1", "TGFBR3", "MAF",
        "IKZF2", "KLF2", "ZEB2", "RGCC", "DUSP2", "DUSP5",
        "NR4A3", "FOS", "BTG2", "RGS2", "RGS16"
    ]
}

gene_to_modules = {}
for module, genes in modules.items():
    for g in genes:
        gene_to_modules.setdefault(g, []).append(module)

def assign_module(g):
    mods = gene_to_modules.get(g, [])
    if len(mods) == 0:
        return "Other"
    return "|".join(mods)

meta["module"] = meta["gene"].apply(assign_module)
meta["rank_log2FC"] = np.arange(1, len(meta) + 1)
meta["minus_log10_padj"] = -np.log10(meta["padj"].replace(0, 1e-300))

# Weighted contribution score for plotting
meta["signature_weight"] = meta["log2FC"] * meta["minus_log10_padj"].clip(upper=300)

# Top subsets
top50 = meta.head(50).copy()
top100 = meta.head(100).copy()
top200 = meta.head(200).copy()

# Module summaries
def module_summary(df, label):
    rows = []
    for module in sorted(df["module"].unique()):
        sub = df[df["module"] == module]
        rows.append({
            "set": label,
            "module": module,
            "genes": len(sub),
            "mean_log2FC": sub["log2FC"].mean(),
            "median_log2FC": sub["log2FC"].median(),
            "mean_score": sub["score"].mean(),
            "mean_minus_log10_padj": sub["minus_log10_padj"].replace(np.inf, 300).mean(),
            "sum_signature_weight": sub["signature_weight"].replace(np.inf, 0).sum(),
            "top_genes": ", ".join(sub["gene"].head(15).tolist())
        })
    return pd.DataFrame(rows)

module_all = pd.concat([
    module_summary(top50, "top50"),
    module_summary(top100, "top100"),
    module_summary(top200, "top200"),
    module_summary(meta, "all_metastasis_up")
], ignore_index=True)

# Gene-level high-confidence table
high_conf = meta[
    (meta["padj"] < 0.05) &
    (meta["log2FC"] > 1.0)
].copy()

# Pathway context if available
pathways = None
if os.path.exists(pathway_file):
    pathways = pd.read_csv(pathway_file)

# Save CSVs
meta.to_csv(os.path.join(FIG7, "7A_metastasis_signature_all_genes_annotated.csv"), index=False)
top50.to_csv(os.path.join(FIG7, "7A_metastasis_signature_top50_annotated.csv"), index=False)
top100.to_csv(os.path.join(FIG7, "7A_metastasis_signature_top100_annotated.csv"), index=False)
module_all.to_csv(os.path.join(FIG7, "7A_module_summary.csv"), index=False)
high_conf.to_csv(os.path.join(FIG7, "7A_high_confidence_metastasis_genes.csv"), index=False)

# Matrix for later visualization: rank x module
module_rank_matrix = pd.crosstab(
    pd.cut(meta["rank_log2FC"], bins=[0, 25, 50, 100, 200, 500], labels=["1-25", "26-50", "51-100", "101-200", "201-500"]),
    meta["module"]
)
module_rank_matrix.to_csv(os.path.join(FIG7, "7A_module_by_rank_bin_matrix.csv"))

# TXT report
report_path = os.path.join(FIG7, "7A_SIGNATURE_NUMERIC_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 7A METASTASIS SIGNATURE NUMERIC REPORT\n\n")

    f.write("1. Input DE file\n")
    f.write(de_file + "\n\n")

    f.write("2. Metastasis-up gene counts\n")
    f.write("Total metastasis-up genes: " + str(len(meta)) + "\n")
    f.write("High-confidence genes log2FC > 1 and padj < 0.05: " + str(len(high_conf)) + "\n\n")

    f.write("3. Top 50 metastasis-up genes with modules\n")
    cols = ["rank_log2FC", "gene", "log2FC", "score", "padj", "minus_log10_padj", "module"]
    f.write(top50[cols].to_string(index=False))
    f.write("\n\n")

    f.write("4. Module summary\n")
    f.write(module_all.to_string(index=False))
    f.write("\n\n")

    f.write("5. Module by rank-bin matrix\n")
    f.write(module_rank_matrix.to_string())
    f.write("\n\n")

    if pathways is not None:
        f.write("6. Top metastasis pathways from Figure 5D\n")
        keep = [c for c in ["Term", "Adjusted P-value", "Combined Score"] if c in pathways.columns]
        f.write(pathways[keep].head(20).to_string(index=False))
        f.write("\n\n")

    f.write("7. Suggested visualization variables\n")
    f.write("x = rank_log2FC\n")
    f.write("y = log2FC or signature_weight\n")
    f.write("color = module\n")
    f.write("size = minus_log10_padj capped at 300\n")
    f.write("labels = top genes per module\n")

print("DONE 7A NUMERIC")
print(report_path)
print(open(report_path, encoding="utf-8").read())