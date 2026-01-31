###-------一,PERMANOVA--------#######
########################
## 0) Paths & load()s
########################
rm(list=ls())

#path_adonis <- 
#out_dir     <- 

load(path_adonis)

########################
## 1) Libraries
########################
library(ggplot2)
library(dplyr)
library(scales)

########################
## 2) Data prep
########################
adonis_df$Pheno <- rownames(adonis_df)

adonis_df <- adonis_df %>%
  filter(
    Pheno != "year",
    Pheno != "Diagnosis2",
    Pheno != "ApoE4"
  )

adonis_df <- adonis_df %>%
  mutate(
    group = case_when(
      Pheno %in% c("ci_per", "gbi_per", "teeth_lost", "teeth_clean",
                   "pd_mean", "pd5_per", "cal_mean", "cal4_per",
                   "periodontitis3") ~ "Periodontal Status",
      Pheno %in% c("MMSE", "MoCA", "Diagnosis") ~ "Cognitive Function",
      Pheno %in% c("diabetes", "Hyperlipidemia", "hypertension") ~ "Medical History",
      Pheno %in% c("age", "sex", "year", "edu", "BMI", "smoke", "alcohol") ~ "Socio-demographic & Lifestyle",
      TRUE ~ NA_character_
    ),
    group = factor(
      group,
      levels = c(
        "Periodontal Status",
        "Cognitive Function",
        "Medical History",
        "Socio-demographic & Lifestyle"
      )
    )
  )

adonis_df <- adonis_df[!is.na(adonis_df$group), ]

adonis_df <- adonis_df %>%
  mutate(
    Pheno = recode_factor(
      Pheno,
      teeth_lost      = "MT",
      gbi_per         = "GBI Sites (%)",
      ci_per          = "CI Sites (%)",
      pd_mean         = "Mean PD",
      cal_mean        = "Mean CAL",
      cal4_per        = "CAL ≥ 4 (%)",
      pd5_per         = "PD ≥ 5 (%)",
      periodontitis3  = "Periodontitis",
      teeth_clean     = "Brushing Frequency",
      Diagnosis       = "Cognitive Status",
      age             = "Age",
      alcohol         = "Alcohol drinking",
      edu             = "Education",
      sex             = "Sex",
      smoke           = "Smoking",
      diabetes        = "Diabetes",
      Hyperlipidemia  = "Hyperlipidemia",
      hypertension    = "Hypertension"
    )
  )

adonis_df <- adonis_df %>%
  mutate(
    R2_percent = R2 * 100,
    anno = case_when(
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

alc_g_color <- c("#C1D09D", "#F8CBC9", "#A4C2B7", "#F29E94")

desired_group_order <- c(
  "Medical History",
  "Cognitive Function",
  "Periodontal Status",
  "Socio-demographic & Lifestyle"
)

adonis_df <- adonis_df %>%
  mutate(group = factor(group, levels = desired_group_order)) %>%
  arrange(group, R2_percent) %>%
  mutate(Pheno = factor(Pheno, levels = Pheno))

########################
## 3) Plot
########################
explain_r <- ggplot(adonis_df, aes(y = Pheno, x = R2)) +
  geom_col(aes(color = group, fill = group), width = 0.8) +
  geom_text(aes(label = anno), size = 6, vjust = 0.5, hjust = -0.2) +
  theme_classic() +
  scale_fill_manual(values = alc_g_color) +
  scale_color_manual(values = alc_g_color) +
  scale_x_continuous(
    limits = c(0, 0.015),
    breaks = c(0, 0.01, 0.015),
    labels = c("0%", "1%", "1.5%"),
    expand = c(0.01, 0)
  ) +
  labs(
    title = "PERMANOVA",
    x = "R² in PERMANOVA, %"
  ) +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13, hjust = 0.6),
    axis.text.y = element_text(size = 13),
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.35)
  ) +
  guides(
    fill  = guide_legend(reverse = TRUE),
    color = guide_legend(reverse = TRUE)
  )

print(explain_r)

########################
## 4) Save
########################
ggsave(
  filename = file.path(out_dir, "adonis2.pdf"),
  plot     = explain_r,
  width    = 8,
  height   = 6,
  device   = cairo_pdf
)

ggsave(
  filename = file.path(out_dir, "adonis2.png"),
  plot     = explain_r,
  width    = 8,
  height   = 6,
  dpi      = 300
)













###----二. Mantel test-----#######


########################
## 0) Paths & load()s
########################
rm(list=ls())

#path_meta <- 
#path_otu  <- 
#path_tax  <- 
#path_mantel <- 

out_dir   <- 
filename  <- 

load(path_meta)
load(path_otu)
load(path_tax)
load(path_mantel)

########################
## 1) Libraries
########################
library(linkET)
library(FD)
library(dplyr)
library(ggplot2)
library(devtools)
library(grid)
library(gridExtra)
library(patchwork)

########################
## 2) Data prep
########################
envs <- env
mantel_times$env <- factor(
  mantel_times$env,
  levels = c("teeth_lost", "gbi_per", "ci_per", "pd_mean", "cal_mean"),
  labels = c("MT", "GBI sites (%)", "CI sites (%)", "Mean PD", "Mean CAL")
)

colnames(envs) <- c(
  "MT", "GBI sites (%)", "Mean PD", "Mean CAL", "CI sites (%)"
)

########################
## 3) Plot Mantel Results
########################
p <- qcorrplot(
  correlate(envs, method = "pearson"),
  type = "lower", diag = FALSE
) +
  geom_square() +
  geom_couple(
    aes(colour = pd, size = rd),
    data = mantel_times,
    curvature = 0.1
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("#3FA9F5", "white", "#e74a32"))(10),
    limits = c(-1, 1), breaks = seq(-1, 1, 0.5)
  ) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(
    size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
    colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
    fill = guide_colorbar(title = "Pearson's r", order = 3)
  ) +
  labs(title = "Mantel test") +
  geom_mark(size = 4, only_mark = TRUE, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, colour = "black") +
  coord_equal(clip = "off") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    plot.margin = margin(6, 40, 6, 6)
  )

print(p)

########################
## 4) Save Mantel plot
########################
png_filename <- file.path(out_dir, paste0(filename, ".png"))
png(png_filename, width = 2000, height = 1350, res = 150)
print(p)
dev.off()

pdf_filename <- file.path(out_dir, paste0(filename, ".pdf"))
pdf(pdf_filename, width = 12, height = 9)
print(p)
dev.off()

########################
## 5) Further Mantel Test Calculation
########################
mantel_times <- mantel_times %>%
  mutate(
    line_type = ifelse(r < 0, "neg", "pos"),
    size_cat  = case_when(
      abs(r) < 0.05 ~ "thin",
      abs(r) < 0.10 ~ "thick",
      TRUE ~ "extra"
    ),
    size_cat = factor(size_cat, levels = c("thin", "thick", "extra"))
  )

qcorrplot(
  correlate(envs, method = "pearson"),
  type = "lower",
  diag = FALSE
) +
  geom_square() +
  geom_couple(
    data = mantel_times,
    aes(linetype = line_type, size = size_cat, colour = pd),
    curvature = 0.1
  ) +
  scale_fill_gradientn(
    name = "Pearson’s r",
    colours = colorRampPalette(c("#3FA9F5", "white", "#e74a32"))(10),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5)
  ) +
  scale_linetype_manual(
    name = "Sign of Mantel’s r",
    values = c(neg = "dashed", pos = "solid")
  ) +
  scale_size_manual(
    name = "Magnitude of Mantel’s r",
    values = c(thin = 0.5, thick = 1, extra = 2),
    guide = guide_legend(override.aes = list(linetype = c("dashed", "solid", "solid")))
  ) +
  scale_colour_manual(
    name = "Mantel’s p",
    values = color_pal(3)
  ) +
  guides(
    linetype = guide_legend(order = 1),
    size = guide_legend(order = 2),
    colour = guide_legend(order = 3),
    fill = guide_colorbar(order = 4)
  ) +
  labs(title = "Mantel test") +
  geom_mark(
    size = 4,
    only_mark = TRUE,
    sig_level = c(0.05, 0.01, 0.001),
    sig_thres = 0.05,
    colour = "black"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
