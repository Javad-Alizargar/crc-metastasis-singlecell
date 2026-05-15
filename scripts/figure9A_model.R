library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# ----------------------------
# Build graph
# ----------------------------
g <- grViz("
digraph CRC_model {

graph [layout = dot, rankdir = TB]

node [shape = box, style = filled, fontname = Helvetica, fontsize=18]

Primary [label = 'Primary CRC\n(epithelial + structural program)', fillcolor = '#FDD0A2']

Stress [label = 'Stress-response activation\n(HSP / protein folding)', fillcolor = '#FB6A4A']

Immune [label = 'Immune remodeling', fillcolor = '#9ECAE1']

Myeloid [label = 'Myeloid / macrophage\nlate immune state', fillcolor = '#6BAED6']

Metastasis [label = 'Metastatic CRC ecosystem', fillcolor = '#CB181D']

Drug [label = 'Drug perturbation landscape\n(mostly mimic stress)', fillcolor = '#DDDDDD']

Target [label = 'Therapeutic target axis\n(HSPA1A, HSPA6, HSPH1, DNAJB1)', fillcolor = '#74C476']

# Main trajectory
Primary -> Stress -> Immune -> Myeloid -> Metastasis

# Branches
Stress -> Drug
Stress -> Target

}
")

# ----------------------------
# SHOW in RStudio viewer
# ----------------------------
print(g)

# ----------------------------
# SAVE high-quality output
# ----------------------------
out_base <- "D:/CRC_META_FULL_SCVI/FIGURE_9_CONCEPTUAL_MODEL/Figure9A_model"

svg <- export_svg(g)
writeLines(svg, paste0(out_base, ".svg"))

rsvg_png(charToRaw(svg), paste0(out_base, ".png"), width = 2400, height = 3200)
rsvg_pdf(charToRaw(svg), paste0(out_base, ".pdf"))

cat("DONE Figure 9A\n")