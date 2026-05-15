# ============================================================
# Figure 7C — Drug-gene chord diagram
# ============================================================

library(circlize)

base_dir <- "D:/CRC_META_FULL_SCVI"
fig7 <- file.path(base_dir, "FIGURE_7_DRUG_REVERSAL_ADVANCED")

edges <- read.csv(file.path(fig7, "7C_chord_edges.csv"), stringsAsFactors = FALSE)
nodes <- read.csv(file.path(fig7, "7C_nodes.csv"), stringsAsFactors = FALSE)

# Shorten drug names for readability
short_drug <- function(x) {
  x <- gsub("LJP006 ", "", x)
  x <- gsub("LJP005 ", "", x)
  x <- gsub("HS578T ", "HS578T_", x)
  x <- gsub("24H-", "24H_", x)
  x <- gsub("3H-", "3H_", x)
  x <- gsub("withaferin-a", "withaferin", x)
  x <- gsub("tanespimycin", "tanesp.", x)
  x <- gsub("geldanamycin", "geldan.", x)
  x <- gsub("puromycin", "puro.", x)
  x <- gsub("menadione", "menad.", x)
  x <- gsub("sanguinarine", "sang.", x)
  return(x)
}

edges$drug_short <- short_drug(edges$drug_term)

# Matrix drug x gene
mat <- xtabs(weight ~ drug_short + gene, data = edges)

# Keep row/column order by total link number
mat <- mat[order(rowSums(mat), decreasing = TRUE), order(colSums(mat), decreasing = TRUE)]

# Colors
module_colors <- c(
  Stress_heatshock = "#E41A1C",
  Signaling_regulatory = "#984EA3",
  Other = "#BDBDBD",
  Metabolic_lipid = "#4DAF4A",
  Antigen_presentation_myeloid = "#A65628",
  Cytotoxic_NK_T = "#377EB8",
  Cytokine_inflammatory = "#FF7F00"
)

gene_module <- unique(edges[, c("gene", "module")])
gene_module <- gene_module[!duplicated(gene_module$gene), ]
gene_module_map <- setNames(gene_module$module, gene_module$gene)

sectors <- c(rownames(mat), colnames(mat))

grid_col <- setNames(rep("#444444", length(sectors)), sectors)

for (g in colnames(mat)) {
  mod <- gene_module_map[[g]]
  if (!is.null(mod) && mod %in% names(module_colors)) {
    grid_col[g] <- module_colors[[mod]]
  } else {
    grid_col[g] <- "#BDBDBD"
  }
}

# Drug sectors neutral dark gray
for (d in rownames(mat)) {
  grid_col[d] <- "#333333"
}

# Link colors by target gene module
link_col <- matrix("#999999", nrow = nrow(mat), ncol = ncol(mat))
rownames(link_col) <- rownames(mat)
colnames(link_col) <- colnames(mat)

for (g in colnames(mat)) {
  mod <- gene_module_map[[g]]
  col <- ifelse(mod %in% names(module_colors), module_colors[[mod]], "#BDBDBD")
  link_col[, g] <- col
}

# Save high-res TIFF
tiff(
  file.path(fig7, "7C_drug_gene_chord.tiff"),
  width = 9000,
  height = 9000,
  res = 700,
  compression = "lzw"
)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(4, nrow(mat) - 1), 18, rep(2, ncol(mat) - 1), 18),
  track.margin = c(0.005, 0.005),
  canvas.xlim = c(-1.25, 1.25),
  canvas.ylim = c(-1.25, 1.25)
)

chordDiagram(
  mat,
  grid.col = grid_col,
  col = link_col,
  transparency = 0.20,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.13),
  link.sort = TRUE,
  link.decreasing = TRUE
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")

    cex_use <- ifelse(sector.name %in% rownames(mat), 0.78, 0.72)

    circos.text(
      mean(xlim),
      ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = cex_use
    )
  },
  bg.border = NA
)

circos.clear()
dev.off()

# Save PDF too
pdf(file.path(fig7, "7C_drug_gene_chord.pdf"), width = 10, height = 10)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(4, nrow(mat) - 1), 18, rep(2, ncol(mat) - 1), 18),
  track.margin = c(0.005, 0.005),
  canvas.xlim = c(-1.25, 1.25),
  canvas.ylim = c(-1.25, 1.25)
)

chordDiagram(
  mat,
  grid.col = grid_col,
  col = link_col,
  transparency = 0.20,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.13),
  link.sort = TRUE,
  link.decreasing = TRUE
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")

    cex_use <- ifelse(sector.name %in% rownames(mat), 0.78, 0.72)

    circos.text(
      mean(xlim),
      ylim[1],
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = cex_use
    )
  },
  bg.border = NA
)

circos.clear()
dev.off()

cat("DONE 7C chord\n")
cat(file.path(fig7, "7C_drug_gene_chord.tiff"), "\n")