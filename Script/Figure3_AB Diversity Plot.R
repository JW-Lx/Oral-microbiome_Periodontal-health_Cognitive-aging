###--------α diversity------###############

########################
## 0) Paths & load()s
########################
rm(list=ls())

#path_dx   <- 
#path_stat <- 

load(path_dx)
load(path_stat)

########################
## 1) Libraries
########################
library(tidyverse)
library(gghalves)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggsignif)
library(ggpubr)
library(showtext)
library(systemfonts)

########################
## 2) Data prep
########################
dx <- dx[!is.na(dx$Group), ]
data <- dx[, c("Group", "Shannon.Wiener")]

df <- data %>%
  mutate(Group = factor(Group, levels = c("Mild", "Moderate", "Severe")))

str(stat.test_Per)
summary(df)

########################
## 3) Plot
########################
font_path <- systemfonts::match_font("Arial")$path
font_add("arial", font_path)
showtext_auto()

theme_set(theme_classic(base_family = "arial"))

ordercolors <- c("#2FA128", "#2D83BE", "#E71D36")

p <- ggplot(df, aes(x = Group, y = Shannon.Wiener, fill = Group)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, width = 0.2, linewidth = 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 3, color = "white") +
  scale_fill_manual(
    values = ordercolors,
    breaks = c("Mild", "Moderate", "Severe"),
    labels = c("None/Mild Periodontitis", "Moderate Periodontitis", "Severe Periodontitis")
  ) +
  scale_y_continuous(limits = c(2.3, 5.5), expand = c(0, 0)) +
  scale_x_discrete(labels = c("None/Mild", "Moderate", "Severe")) +
  labs(y = "Shannon    Index", x = "Periodontitis    Group") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16, color = "black"),
    axis.text  = element_text(size = 13, color = "black")
  )

print(p)

stat_shannon <- stat.test_Per %>%
  filter(Index == "Shannon.Wiener") %>%
  mutate(
    y.position = c(5.3, 5.4, 5.55),
    label.y    = c(5.3, 5.4, 5.55)
  )

p1 <- p +
  stat_pvalue_manual(
    data = stat_shannon,
    mapping = aes(
      x = group1,
      xend = group2,
      y = y.position,
      label = label
    ),
    inherit.aes = FALSE,
    tip.length = 0.015,
    bracket.size = 0.5,
    bracket.shorten = 0.01,
    hide.ns = TRUE,
    size = 6
  )

print(p1)

########################
## 4) Outliers summary
########################
df %>%
  group_by(Group) %>%
  summarise(
    Q1        = quantile(Shannon.Wiener, .25),
    Q3        = quantile(Shannon.Wiener, .75),
    IQR       = IQR(Shannon.Wiener),
    lower     = Q1 - 1.5 * IQR,
    upper     = Q3 + 1.5 * IQR,
    n_outlier = sum(Shannon.Wiener < lower | Shannon.Wiener > upper),
    total     = n(),
    prop      = n_outlier / total
  )



###--------β diversity------###############


########################
## 0) Paths & load()s
########################
rm(list=ls())

#path_meta <- 
#path_otu  <- 
#path_tax  <- 

#out_dir   <- 
#filename  <- 

load(path_meta)
load(path_otu)
load(path_tax)

########################
## 1) Libraries
########################
library(vegan)
library(ape)
library(ggplot2)
library(phyloseq)
library(gplots)
library(RColorBrewer)
library(ggExtra)


########################
## 2) PCoA plot
########################
pcoa <- ordinate(ps_norm, method = "PCoA", distance = "bray")
var_exp <- round(100 * pcoa$values$Relative_eig[1:2], 1)
x_lab <- paste0("PcoA1 [", var_exp[1], "%]")
y_lab <- paste0("PcoA2 [", var_exp[2], "%]")

p <- plot_ordination(ps_norm, pcoa, color = "Diagnosis") +
  stat_ellipse(geom = "polygon", alpha = 0.1, aes(fill = Diagnosis)) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Cognitive Status",
    values = c("CN" = "#2FA128", "MCI" = "#2D83BE", "D" = "#E71D36")
  ) +
  scale_fill_manual(
    name = "Cognitive Status",
    values = c("CN" = "#2FA128", "MCI" = "#2D83BE", "D" = "#E71D36")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = x_lab, y = y_lab)

convert_p_to_stars <- function(p) {
  if (!is.na(p) && p < 0.01) {
    "**"
  } else if (!is.na(p) && p < 0.05) {
    "*"
  } else {
    ""
  }
}

sig_text <- c(
  sprintf("CN vs MCI: R² = %.3f%s",
          result_df["CN_vs_MCI", "R2"],
          convert_p_to_stars(result_df["CN_vs_MCI", "p.adj"])),
  sprintf("CN vs D:   R² = %.3f%s",
          result_df["CN_vs_D", "R2"],
          convert_p_to_stars(result_df["CN_vs_D", "p.adj"])),
  sprintf("MCI vs D:  R² = %.3f%s",
          result_df["MCI_vs_D", "R2"],
          convert_p_to_stars(result_df["MCI_vs_D", "p.adj"]))
)

p2 <- p +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste(sig_text, collapse = "\n"),
    hjust = -0.1, vjust = 1.2,
    size = 3.8,
    color = "black",
    fontface = "italic"
  )

p_marginal <- ggMarginal(p2, type = "density", groupColour = TRUE, groupFill = TRUE)
print(p_marginal)

########################
## 5) Save
########################
png_filename <- file.path(out_dir, paste0(filename, ".png"))
png(png_filename, width = 2000, height = 1350, res = 150)
print(p_marginal)
dev.off()

pdf_filename <- file.path(out_dir, paste0(filename, ".pdf"))
pdf(pdf_filename, width = 12, height = 9)
print(p_marginal)
dev.off()