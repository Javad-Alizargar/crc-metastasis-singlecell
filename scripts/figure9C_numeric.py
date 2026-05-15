# -*- coding: utf-8 -*-

import os
import pandas as pd

BASE_IMMUNE = r"D:\CRC_META_FULL_SCVI\FIGURE_6_TRAJECTORY_IMMUNE"
FIG9 = r"D:\CRC_META_FULL_SCVI\FIGURE_9_CONCEPTUAL_MODEL"
os.makedirs(FIG9, exist_ok=True)

quartile_file = os.path.join(BASE_IMMUNE, "STAT_pseudotime_quartile_celltype_percent.csv")
celltype_summary_file = os.path.join(BASE_IMMUNE, "STAT_celltype_pseudotime_summary.csv")
within_ct_file = os.path.join(BASE_IMMUNE, "STAT_within_celltype_primary_vs_metastasis_tests.csv")

quartile = pd.read_csv(quartile_file)
celltype_summary = pd.read_csv(celltype_summary_file)
within_ct = pd.read_csv(within_ct_file)

# ------------------------------------------------------------
# Long format: quartile -> cell type percent
# ------------------------------------------------------------
long = quartile.melt(
    id_vars=["pseudotime_quartile"],
    var_name="cell_type",
    value_name="percent"
)

# ------------------------------------------------------------
# Ecosystem state labels
# ------------------------------------------------------------
def ecosystem_state(row):
    q = str(row["pseudotime_quartile"])
    ct = str(row["cell_type"])

    if q == "Q1_early":
        return "Early immune state"
    if q in ["Q2", "Q3"]:
        return "Intermediate immune state"
    if q == "Q4_late":
        return "Late myeloid/macrophage state"
    return "Unassigned"

long["ecosystem_state"] = long.apply(ecosystem_state, axis=1)

# ------------------------------------------------------------
# Add broad immune class
# ------------------------------------------------------------
broad_map = {
    "Plasma_cells": "Adaptive humoral",
    "Cytotoxic_T_NK": "Cytotoxic lymphoid",
    "Macrophage": "Innate myeloid",
    "Myeloid": "Innate myeloid"
}
long["broad_class"] = long["cell_type"].map(broad_map).fillna("Other")

# ------------------------------------------------------------
# Add mean pseudotime per cell type
# ------------------------------------------------------------
ct_mean = dict(zip(celltype_summary["cell_type"], celltype_summary["mean"]))
long["celltype_mean_pseudotime"] = long["cell_type"].map(ct_mean)

# ------------------------------------------------------------
# Add metastasis-primary shift per cell type
# ------------------------------------------------------------
shift_map = dict(
    zip(
        within_ct["cell_type"],
        within_ct["mean_difference_metastasis_minus_primary"]
    )
)
long["metastasis_minus_primary_pseudotime_shift"] = long["cell_type"].map(shift_map)

# ------------------------------------------------------------
# Create simplified alluvial table
# ------------------------------------------------------------
alluvial = long.copy()
alluvial["from_stage"] = alluvial["pseudotime_quartile"]
alluvial["middle_class"] = alluvial["broad_class"]
alluvial["to_cell_type"] = alluvial["cell_type"]
alluvial["weight"] = alluvial["percent"]

alluvial = alluvial[
    [
        "from_stage",
        "middle_class",
        "to_cell_type",
        "ecosystem_state",
        "weight",
        "celltype_mean_pseudotime",
        "metastasis_minus_primary_pseudotime_shift"
    ]
].copy()

# ------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------
state_summary = (
    long.groupby(["pseudotime_quartile", "ecosystem_state", "broad_class"], observed=True)["percent"]
    .sum()
    .reset_index()
    .sort_values(["pseudotime_quartile", "percent"], ascending=[True, False])
)

celltype_by_quartile = long.pivot_table(
    index="pseudotime_quartile",
    columns="cell_type",
    values="percent",
    aggfunc="sum"
).reset_index()

# ------------------------------------------------------------
# Save files
# ------------------------------------------------------------
alluvial.to_csv(os.path.join(FIG9, "Figure9C_alluvial_input.csv"), index=False)
state_summary.to_csv(os.path.join(FIG9, "Figure9C_state_summary.csv"), index=False)
celltype_by_quartile.to_csv(os.path.join(FIG9, "Figure9C_celltype_by_quartile.csv"), index=False)

# ------------------------------------------------------------
# TXT report
# ------------------------------------------------------------
report = os.path.join(FIG9, "Figure9C_NUMERIC_REPORT.txt")

with open(report, "w", encoding="utf-8") as f:
    f.write("FIGURE 9C IMMUNE ECOSYSTEM TRANSITION NUMERIC REPORT\n\n")

    f.write("1. Input files\n")
    f.write(quartile_file + "\n")
    f.write(celltype_summary_file + "\n")
    f.write(within_ct_file + "\n\n")

    f.write("2. Cell type percent by pseudotime quartile\n")
    f.write(celltype_by_quartile.to_string(index=False))
    f.write("\n\n")

    f.write("3. Ecosystem state summary\n")
    f.write(state_summary.to_string(index=False))
    f.write("\n\n")

    f.write("4. Cell type mean pseudotime\n")
    f.write(celltype_summary.to_string(index=False))
    f.write("\n\n")

    f.write("5. Within-cell-type metastasis minus primary shift\n")
    f.write(within_ct.to_string(index=False))
    f.write("\n\n")

    f.write("6. Key interpretation\n")
    f.write("Q1_early is plasma-cell enriched.\n")
    f.write("Q2 and Q3 are cytotoxic T/NK enriched intermediate states.\n")
    f.write("Q4_late is macrophage/myeloid enriched.\n")
    f.write("This supports an immune ecosystem transition from adaptive/humoral and cytotoxic states toward innate myeloid/macrophage states.\n\n")

    f.write("7. Files generated\n")
    for x in sorted(os.listdir(FIG9)):
        f.write(x + "\n")

print("DONE FIGURE 9C NUMERIC")
print(report)
print(open(report, encoding="utf-8").read())