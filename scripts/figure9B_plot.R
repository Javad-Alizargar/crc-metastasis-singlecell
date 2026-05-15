library(tidyverse)

# Load data
df <- read.csv("D:/CRC_META_FULL_SCVI/FIGURE_9_CONCEPTUAL_MODEL/Figure9B_evidence_matrix_long.csv")

# Order factors (IMPORTANT for clean figure)
df$short_label <- factor(df$short_label, levels = c(
  "Stress axis",
  "Immune activation",
  "Myeloid late state",
  "Plasma early state",
  "Primary structure",
  "Tumor stability",
  "Drug mimicry",
  "Target axis"
))

df$Evidence_layer_label <- factor(df$Evidence_layer_label, levels = c(
  "Differential expression",
  "Pathway enrichment",
  "Score overlays",
  "Pseudotime",
  "Perturbation network"
))

df$Strength_label <- factor(df$Strength_label, levels = c(
  "None", "Weak", "Moderate", "Strong"
))

# Color palette (VERY IMPORTANT)
colors <- c(
  "None" = "#F5F5F5",
  "Weak" = "#BDBDBD",
  "Moderate" = "#6BAED6",
  "Strong" = "#E41A1C"
)

# Plot
p <- ggplot(df, aes(x = Evidence_layer_label, y = short_label, fill = Strength_label)) +
  geom_tile(color = "white", size = 1.2) +

  # Add numbers inside tiles (clean!)
  geom_text(aes(label = Strength), size = 5, fontface = "bold") +

  scale_fill_manual(values = colors) +

  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +

  labs(
    x = "",
    y = "",
    fill = "Evidence strength"
  )

# Save
ggsave(
  "D:/CRC_META_FULL_SCVI/FIGURE_9_CONCEPTUAL_MODEL/Figure9B_evidence_heatmap.tiff",
  p,
  width = 12,
  height = 7,
  dpi = 600
)

ggsave(
  "D:/CRC_META_FULL_SCVI/FIGURE_9_CONCEPTUAL_MODEL/Figure9B_evidence_heatmap.pdf",
  p,
  width = 12,
  height = 7
)

print(p)