########################
## 0) Paths & load()s
########################
rm(list=ls())

#work_dir  <- 
#data_xlsx <- 
#pdf_file  <- 

setwd(work_dir)
data <- read_xlsx(data_xlsx)

########################
## 1) Libraries
########################
library(readxl)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(eulerr)
library(grid)
library(dplyr)
library(stringr)

########################
## 2) Data prep
########################
data <- data %>%
  mutate(
    ID = if_else(Group == "Genus", str_remove(feature, "^g_"), ID),
    ID = factor(ID, levels = ID)
  )

colnames(data) <- c(
  "Group",
  "feature",
  "Mean CAL",
  "CAL ≥ 4 %",
  "CI (%)",
  "PD ≥ 5 %",
  "ModP",
  "SevP",
  "Mean CAL_qval",
  "CAL ≥ 4 %_qval",
  "CI (%)_qval",
  "PD ≥ 5 %_qval",
  "ModP_qval",
  "SevP_qval",
  "ID",
  "abundance"
)

data_normalized <- data %>%
  group_by(Group) %>%
  arrange(abundance, .by_group = TRUE) %>%
  ungroup()

data_matrix <- data_normalized %>%
  select("Mean CAL", "CAL ≥ 4 %", "CI (%)", "PD ≥ 5 %", "ModP", "SevP") %>%
  as.matrix()

data_matrix_heat <- data_matrix[, ncol(data_matrix):1]
rownames(data_matrix_heat) <- data_normalized$ID
rownames(data_matrix) <- data_normalized$ID

dat_matrix_rev <- as.data.frame(data_matrix)
dt <- as.data.frame(data_normalized)
rownames(dt) <- dt$ID
all.equal(rownames(dat_matrix_rev), rownames(dt))

col_map <- c(
  `Mean CAL`    = "Mean CAL_qval",
  `CAL ≥ 4 %`   = "CAL ≥ 4 %_qval",
  `CI (%)`      = "CI (%)_qval",
  `PD ≥ 5 %`    = "PD ≥ 5 %_qval",
  ModP          = "ModP_qval",
  SevP          = "SevP_qval"
)

sig_mat <- dat_matrix_rev
sig_mat[] <- lapply(sig_mat, as.character)

for (heat_col in names(col_map)) {
  qcol <- col_map[[heat_col]]
  qvec <- dt[rownames(sig_mat), qcol]
  sig_mat[[heat_col]] <- ifelse(
    is.na(qvec), "",
    ifelse(qvec < 0.01, "**",
           ifelse(qvec < 0.05, "*",
                  ifelse(qvec < 0.10, "+", "")))
  )
}

rng <- quantile(as.matrix(data[, 3:8]), probs = c(0.1, 0.9), na.rm = TRUE)

green_pink <- colorRamp2(
  c(-0.15, 0, 0.4),
  c("#b2182b", "#e4e7ec", "#2166ac")
)

group_colors <- c("Genus" = "#72A6C6", "Pathway" = "#CA8985")

red_black <- colorRamp2(
  seq(0, 1, length.out = 6),
  rev(c("#fff1e7", "#f99c70", "#dd4139", "#97044d", "#4f2043", "black"))
)

data_normalized$Group <- factor(data_normalized$Group, levels = c("Genus", "Pathway"))

########################
## 3) Plot
########################
par(family = "Arial")
pdf(pdf_file, width = 15, height = 15)

circos.clear()
circos.par(
  start.degree = 60,
  gap.after = c(rep(5, length(levels(data_normalized$Group)) - 1), 30),
  track.margin = c(0, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

font_vec <- ifelse(data_normalized$Group == "Genus", 3, 1)

circos.heatmap(
  data_matrix_heat,
  split = data_normalized$Group,
  cluster = FALSE,
  bg.border = "black",
  bg.lwd = 0.2,
  cell.border = "white",
  cell.lwd = 0.2,
  rownames.side = "outside",
  rownames.cex = 0.6,
  rownames.font = font_vec,
  col = green_pink,
  track.height = 0.25
)

for (sector.index in get.all.sector.index()) {
  idx <- which(data_normalized$Group == sector.index)
  values <- sig_mat[idx, , drop = FALSE]
  xlim <- get.cell.meta.data("xlim", sector.index)
  ylim <- get.cell.meta.data("ylim", sector.index)
  n_row <- nrow(values)
  n_col <- ncol(values)
  
  x_points <- seq(xlim[1], xlim[2], length.out = n_row + 1)
  y_points <- seq(ylim[1], ylim[2], length.out = n_col + 1)
  
  for (i in 1:n_row) {
    x <- (x_points[i] + x_points[i + 1]) / 2
    for (j in 1:n_col) {
      y <- (y_points[j] + y_points[j + 1]) / 2
      circos.text(
        x, y,
        labels = values[i, j],
        sector.index = sector.index,
        facing = "inside",
        adj = c(0.5, 0.5),
        cex = 0.8,
        niceFacing = TRUE
      )
    }
  }
}

circos.track(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
      cn <- colnames(dat_matrix_rev)
      n <- length(cn)
      cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      y_coords <- seq(
        CELL_META$cell.ylim[1] + cell_height / 2,
        CELL_META$cell.ylim[2] - cell_height / 2,
        length.out = n
      )
      
      for (i in 1:n) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
          c(y_coords[i], y_coords[i]),
          col = "black",
          lwd = 2
        )
      }
      
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"),
        y_coords,
        cn,
        cex = 1,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)

circos.track(
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    sector_data <- data_normalized[data_normalized$Group == CELL_META$sector.index, ]
    
    for (i in 1:nrow(sector_data)) {
      circos.points(
        CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data),
        0.5,
        pch = 18,
        cex = 1.8,
        col = red_black(sector_data$abundance[i])
      )
    }
    
    if (CELL_META$sector.numeric.index == length(levels(data_normalized$Group))) {
      circos.lines(
        c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
        c(0.5, 0.5),
        col = "black",
        lwd = 2
      )
      circos.text(
        CELL_META$cell.xlim[2] + convert_x(1.5, "mm"),
        0.5,
        "Abundance",
        cex = 1,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)

circos.track(
  ylim = c(0, 1),
  track.height = 0.065,
  bg.col = adjustcolor(group_colors[levels(data_normalized$Group)], alpha.f = 0.3),
  bg.border = "black",
  bg.lwd = 0.2,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.75,
      CELL_META$sector.index,
      facing = "bending.inside",
      cex = 1.5,
      adj = c(0.5, 0)
    )
  }
)

heatmap_legend <- Legend(
  title = "Coef",
  col_fun = green_pink,
  at = c(-0.15, 0, 0.40),
  labels = c("-0.15", "0", "0.40"),
  title_position = "leftcenter-rot",
  title_gp = gpar(fontsize = 14),
  labels_gp = gpar(fontsize = 14)
)

qval_legend <- Legend(
  title = "abundance",
  col_fun = red_black,
  at = c(1, 0),
  labels = c("Low", "High"),
  title_position = "leftcenter-rot",
  title_gp = gpar(fontsize = 14),
  labels_gp = gpar(fontsize = 14)
)

draw(
  heatmap_legend,
  x = unit(0.95, "npc") - unit(5, "mm"),
  y = unit(0.9, "npc") - unit(5, "mm"),
  just = c("right", "top")
)

draw(
  qval_legend,
  x = unit(0.95, "npc") - unit(5, "mm"),
  y = unit(0.20, "npc") - unit(5, "mm"),
  just = c("right", "top")
)

n_genus <- data %>% filter(Group == "Genus") %>% nrow()
n_pathway <- data %>% filter(Group == "Pathway") %>% nrow()

pushViewport(viewport(x = 0.5, y = 0.5, width = 1, height = 1, just = c("centre", "centre")))

r0 <- 0.08
x_left <- 0.42
x_right <- 0.58

grid.circle(x = x_left, y = 0.52, r = r0, gp = gpar(fill = group_colors["Pathway"], col = NA, alpha = 0.6))
grid.text(sprintf("%d", n_pathway), x = x_left, y = 0.52, gp = gpar(fontsize = 20, col = "black", fontface = "bold"))

grid.circle(x = x_right, y = 0.52, r = r0, gp = gpar(fill = group_colors["Genus"], col = NA, alpha = 0.6))
grid.text(sprintf("%d", n_genus), x = x_right, y = 0.52, gp = gpar(fontsize = 20, col = "black", fontface = "bold"))

grid.text("Oral Phenotype Associations", x = 0.5, y = 0.40, gp = gpar(fontsize = 20, fontface = "italic"))

upViewport()

dev.off()
