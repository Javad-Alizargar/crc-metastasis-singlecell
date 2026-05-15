library(ggplot2)
library(ggalluvial)
library(dplyr)

fig9_dir <- "D:/CRC_META_FULL_SCVI/FIGURE_9_CONCEPTUAL_MODEL"

df <- read.csv(file.path(fig9_dir, "Figure9C_alluvial_input.csv"))

df$from_stage <- factor(
  df$from_stage,
  levels = c("Q1_early", "Q2", "Q3", "Q4_late")
)

df$middle_class <- factor(
  df$middle_class,
  levels = c("Adaptive humoral", "Cytotoxic lymphoid", "Innate myeloid")
)

df$to_cell_type <- factor(
  df$to_cell_type,
  levels = c("Plasma_cells", "Cytotoxic_T_NK", "Macrophage", "Myeloid")
)

class_cols <- c(
  "Adaptive humoral" = "#984EA3",
  "Cytotoxic lymphoid" = "#4575B4",
  "Innate myeloid" = "#D73027"
)

p <- ggplot(
  df,
  aes(
    axis1 = from_stage,
    axis2 = middle_class,
    axis3 = to_cell_type,
    y = weight
  )
) +
  geom_alluvium(
    aes(fill = middle_class),
    width = 0.16,
    alpha = 0.78,
    knot.pos = 0.45
  ) +
  geom_stratum(
    width = 0.22,
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 4.2,
    fontface = "bold"
  ) +
  scale_x_discrete(
    limits = c("Pseudotime quartile", "Immune class", "Cell type"),
    expand = c(0.08, 0.08)
  ) +
  scale_fill_manual(values = class_cols) +
  labs(
    x = "",
    y = "Percent of cells",
    fill = "Immune class"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(
  file.path(fig9_dir, "Figure9C_immune_ecosystem_alluvial.tiff"),
  p,
  width = 12,
  height = 7,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  file.path(fig9_dir, "Figure9C_immune_ecosystem_alluvial.pdf"),
  p,
  width = 12,
  height = 7,
  bg = "white"
)

print(p)
cat("DONE Figure 9C alluvial\n")