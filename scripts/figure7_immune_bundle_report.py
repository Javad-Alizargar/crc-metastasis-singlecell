# -*- coding: utf-8 -*-

import os
import pandas as pd

BASE = r"D:\CRC_META_FULL_SCVI\FIGURE_6_TRAJECTORY_IMMUNE"
out_file = os.path.join(BASE, "FIGURE7_IMMUNE_FULL_REPORT.txt")

def safe_read_csv(path):
    if os.path.exists(path):
        return pd.read_csv(path)
    return None

report_txt = os.path.join(BASE, "Figure6_IMMUNE_STATISTICAL_REPORT.txt")

label_test = safe_read_csv(os.path.join(BASE, "STAT_label_pseudotime_test.csv"))
celltype_summary = safe_read_csv(os.path.join(BASE, "STAT_celltype_pseudotime_summary.csv"))
score_corr = safe_read_csv(os.path.join(BASE, "STAT_score_pseudotime_correlations.csv"))
quartile = safe_read_csv(os.path.join(BASE, "STAT_pseudotime_quartile_celltype_percent.csv"))
within_ct = safe_read_csv(os.path.join(BASE, "STAT_within_celltype_primary_vs_metastasis_tests.csv"))

with open(out_file, "w", encoding="utf-8") as f:

    f.write("FIGURE 7 IMMUNE TRAJECTORY FULL REPORT\n\n")

    if os.path.exists(report_txt):
        f.write("MAIN STATISTICAL REPORT\n")
        with open(report_txt, encoding="utf-8") as r:
            f.write(r.read())
        f.write("\n\n")

    if label_test is not None:
        f.write("PRIMARY VS METASTASIS PSEUDOTIME\n")
        f.write(label_test.to_string(index=False))
        f.write("\n\n")

    if celltype_summary is not None:
        f.write("CELL TYPE PSEUDOTIME SUMMARY\n")
        f.write(celltype_summary.to_string(index=False))
        f.write("\n\n")

    if score_corr is not None:
        f.write("SCORE PSEUDOTIME CORRELATIONS\n")
        f.write(score_corr.to_string(index=False))
        f.write("\n\n")

    if quartile is not None:
        f.write("CELL TYPE PERCENT BY PSEUDOTIME QUARTILE\n")
        f.write(quartile.to_string(index=False))
        f.write("\n\n")

    if within_ct is not None:
        f.write("WITHIN CELL TYPE PRIMARY VS METASTASIS\n")
        f.write(within_ct.to_string(index=False))
        f.write("\n\n")

    f.write("AUTO SUMMARY\n")

    if label_test is not None:
        row = label_test.iloc[0]
        f.write("Mean difference metastasis minus primary: " + str(row.get("mean_difference_metastasis_minus_primary", "NA")) + "\n")
        f.write("P-value: " + str(row.get("p_value", "NA")) + "\n")
        f.write("Cliffs delta: " + str(row.get("cliffs_delta_metastasis_vs_primary", "NA")) + "\n\n")

    if celltype_summary is not None:
        top_ct = celltype_summary.sort_values("mean", ascending=False).iloc[0]
        f.write("Highest pseudotime cell type: " + str(top_ct["cell_type"]) + " (mean=" + str(top_ct["mean"]) + ")\n\n")

    if score_corr is not None:
        for _, r in score_corr.iterrows():
            f.write(str(r["score"]) + " correlation: " + str(r["spearman_rho_with_pseudotime"]) + "\n")

print("DONE")
print(out_file)