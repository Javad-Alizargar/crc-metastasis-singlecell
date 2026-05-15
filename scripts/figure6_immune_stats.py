import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, kruskal, spearmanr

OUT = r"D:\CRC_META_FULL_SCVI"
FIG6 = os.path.join(OUT, "FIGURE_6_TRAJECTORY_IMMUNE")

adata_path = os.path.join(FIG6, "Figure6_immune_trajectory.h5ad")
adata = sc.read_h5ad(adata_path)

df = adata.obs.copy()
df["label"] = df["label"].astype(str)
df["cell_type"] = df["cell_type"].astype(str)
df["site"] = df["site"].astype(str)

score_cols = [
    "Metastasis_immune_score",
    "Metastasis_stress_score",
    "Primary_structural_score"
]
score_cols = [x for x in score_cols if x in df.columns]

# ------------------------------------------------------------
# Helper
# ------------------------------------------------------------
def cliffs_delta(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    nx = len(x)
    ny = len(y)
    if nx == 0 or ny == 0:
        return np.nan
    # Fast approximate using ranks
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank().values
    rx = ranks[:nx].sum()
    u = rx - nx * (nx + 1) / 2
    delta = (2 * u) / (nx * ny) - 1
    return delta

# ------------------------------------------------------------
# 1. Pseudotime primary vs metastasis
# ------------------------------------------------------------
primary_pt = df.loc[df["label"] == "primary", "dpt_pseudotime"]
meta_pt = df.loc[df["label"] == "metastasis", "dpt_pseudotime"]

u_stat, p_val = mannwhitneyu(primary_pt, meta_pt, alternative="two-sided")
delta = cliffs_delta(meta_pt, primary_pt)

label_test = pd.DataFrame([{
    "comparison": "metastasis_vs_primary",
    "primary_n": len(primary_pt),
    "metastasis_n": len(meta_pt),
    "primary_mean": primary_pt.mean(),
    "metastasis_mean": meta_pt.mean(),
    "primary_median": primary_pt.median(),
    "metastasis_median": meta_pt.median(),
    "mean_difference_metastasis_minus_primary": meta_pt.mean() - primary_pt.mean(),
    "median_difference_metastasis_minus_primary": meta_pt.median() - primary_pt.median(),
    "mannwhitney_u": u_stat,
    "p_value": p_val,
    "cliffs_delta_metastasis_vs_primary": delta
}])

# ------------------------------------------------------------
# 2. Pseudotime by cell type
# ------------------------------------------------------------
celltype_summary = (
    df.groupby("cell_type", observed=True)["dpt_pseudotime"]
    .agg(["count", "mean", "median", "std", "min", "max"])
    .reset_index()
    .sort_values("mean")
)

groups = [g["dpt_pseudotime"].values for _, g in df.groupby("cell_type", observed=True)]
kw_stat, kw_p = kruskal(*groups)

celltype_kw = pd.DataFrame([{
    "test": "kruskal_wallis_cell_type_pseudotime",
    "statistic": kw_stat,
    "p_value": kw_p
}])

# Pairwise cell type tests
pair_rows = []
cell_types = sorted(df["cell_type"].unique())

for i in range(len(cell_types)):
    for j in range(i + 1, len(cell_types)):
        a = cell_types[i]
        b = cell_types[j]
        xa = df.loc[df["cell_type"] == a, "dpt_pseudotime"]
        xb = df.loc[df["cell_type"] == b, "dpt_pseudotime"]
        u, p = mannwhitneyu(xa, xb, alternative="two-sided")
        pair_rows.append({
            "cell_type_A": a,
            "cell_type_B": b,
            "A_mean": xa.mean(),
            "B_mean": xb.mean(),
            "B_minus_A_mean": xb.mean() - xa.mean(),
            "p_value": p,
            "cliffs_delta_B_vs_A": cliffs_delta(xb, xa)
        })

pairwise_celltype = pd.DataFrame(pair_rows)

# ------------------------------------------------------------
# 3. Pseudotime by label and cell type
# ------------------------------------------------------------
label_celltype_summary = (
    df.groupby(["label", "cell_type"], observed=True)["dpt_pseudotime"]
    .agg(["count", "mean", "median", "std"])
    .reset_index()
)

within_ct_rows = []

for ct in cell_types:
    sub = df[df["cell_type"] == ct]
    if set(["primary", "metastasis"]).issubset(set(sub["label"].unique())):
        x = sub.loc[sub["label"] == "primary", "dpt_pseudotime"]
        y = sub.loc[sub["label"] == "metastasis", "dpt_pseudotime"]
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        within_ct_rows.append({
            "cell_type": ct,
            "primary_n": len(x),
            "metastasis_n": len(y),
            "primary_mean": x.mean(),
            "metastasis_mean": y.mean(),
            "mean_difference_metastasis_minus_primary": y.mean() - x.mean(),
            "primary_median": x.median(),
            "metastasis_median": y.median(),
            "median_difference_metastasis_minus_primary": y.median() - x.median(),
            "p_value": p,
            "cliffs_delta_metastasis_vs_primary": cliffs_delta(y, x)
        })

within_celltype_label_tests = pd.DataFrame(within_ct_rows)

# ------------------------------------------------------------
# 4. Correlation of pathway scores with pseudotime
# ------------------------------------------------------------
corr_rows = []

for score in score_cols:
    rho, p = spearmanr(df["dpt_pseudotime"], df[score])
    corr_rows.append({
        "score": score,
        "spearman_rho_with_pseudotime": rho,
        "p_value": p
    })

score_correlations = pd.DataFrame(corr_rows)

# ------------------------------------------------------------
# 5. Score means by pseudotime quartile
# ------------------------------------------------------------
df["pseudotime_quartile"] = pd.qcut(
    df["dpt_pseudotime"],
    q=4,
    labels=["Q1_early", "Q2", "Q3", "Q4_late"]
)

quartile_summary = (
    df.groupby("pseudotime_quartile", observed=True)[["dpt_pseudotime"] + score_cols]
    .mean()
    .reset_index()
)

quartile_celltype = pd.crosstab(
    df["pseudotime_quartile"],
    df["cell_type"],
    normalize="index"
) * 100

# ------------------------------------------------------------
# Save CSVs
# ------------------------------------------------------------
label_test.to_csv(os.path.join(FIG6, "STAT_label_pseudotime_test.csv"), index=False)
celltype_summary.to_csv(os.path.join(FIG6, "STAT_celltype_pseudotime_summary.csv"), index=False)
celltype_kw.to_csv(os.path.join(FIG6, "STAT_celltype_kruskal_test.csv"), index=False)
pairwise_celltype.to_csv(os.path.join(FIG6, "STAT_pairwise_celltype_tests.csv"), index=False)
label_celltype_summary.to_csv(os.path.join(FIG6, "STAT_label_celltype_pseudotime_summary.csv"), index=False)
within_celltype_label_tests.to_csv(os.path.join(FIG6, "STAT_within_celltype_primary_vs_metastasis_tests.csv"), index=False)
score_correlations.to_csv(os.path.join(FIG6, "STAT_score_pseudotime_correlations.csv"), index=False)
quartile_summary.to_csv(os.path.join(FIG6, "STAT_pseudotime_quartile_score_summary.csv"), index=False)
quartile_celltype.to_csv(os.path.join(FIG6, "STAT_pseudotime_quartile_celltype_percent.csv"))

# ------------------------------------------------------------
# TXT report
# ------------------------------------------------------------
report_path = os.path.join(FIG6, "Figure6_IMMUNE_STATISTICAL_REPORT.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("FIGURE 6 IMMUNE TRAJECTORY STATISTICAL REPORT\n\n")

    f.write("1. Primary vs metastasis pseudotime test\n")
    f.write(label_test.to_string(index=False))
    f.write("\n\n")

    f.write("2. Pseudotime by cell type summary\n")
    f.write(celltype_summary.to_string(index=False))
    f.write("\n\n")

    f.write("3. Kruskal-Wallis test across cell types\n")
    f.write(celltype_kw.to_string(index=False))
    f.write("\n\n")

    f.write("4. Within-cell-type primary vs metastasis tests\n")
    f.write(within_celltype_label_tests.to_string(index=False))
    f.write("\n\n")

    f.write("5. Pseudotime score correlations\n")
    f.write(score_correlations.to_string(index=False))
    f.write("\n\n")

    f.write("6. Score means by pseudotime quartile\n")
    f.write(quartile_summary.to_string(index=False))
    f.write("\n\n")

    f.write("7. Cell type percent by pseudotime quartile\n")
    f.write(quartile_celltype.round(2).to_string())
    f.write("\n\n")

print("DONE IMMUNE STATS")
print(report_path)
print(open(report_path, encoding="utf-8").read())