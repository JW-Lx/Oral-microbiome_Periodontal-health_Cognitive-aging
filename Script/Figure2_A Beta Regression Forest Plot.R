########################
## 0) Paths & load()s
########################
rm(list=ls())

#path_dt   <- ""
#path_meta <- ""
#out_dir   <- ""

dt <- readxl::read_excel(path_dt)
load(path_meta)

########################
## 1) Libraries
########################
library(readxl)
library(dplyr)
library(grid)
library(forestploter)
library(stringr)
library(showtext)

########################
## 2) Data prep
########################
df <- metadata[!is.na(metadata$MMSE), ]

dt$`Sample Size` <- c(
  1132, 1114, 1113, 1115, 1115, 1115, 1115, "668/1119", "419/1119",
  1132, 1114, 1113, 1115, 1115, 1115, 1115, "668/1119", "419/1119"
)

dt$Phenotype <- ifelse(is.na(dt$Result), dt$Phenotype, paste0("   ", dt$Phenotype))
dt$`Sample Size` <- ifelse(is.na(dt$Result), dt$`Sample Size`, paste0("      ", dt$`Sample Size`))

dt$Result <- ifelse(is.na(dt$Result), "", dt$Result)
dt$`Beta(95%CI)` <- ifelse(is.na(dt$`Beta(95%CI)`), "", dt$`Beta(95%CI)`)
dt$`Sample Size` <- ifelse(is.na(dt$`Sample Size`), "", dt$`Sample Size`)
dt$se <- (dt$CI_up - dt$CI_low) / (2 * 1.96)

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt1 <- dt[c(1:9), ]
dt2 <- dt[c(10:18), ]

colnames(dt1)[colnames(dt1) == "Estimate"] <- "est_moca"
colnames(dt1)[colnames(dt1) == "CI_low"] <- "low_moca"
colnames(dt1)[colnames(dt1) == "CI_up"] <- "hi_moca"
colnames(dt1)[colnames(dt1) == "Beta(95%CI)"] <- "Beta(95%CI)_moca"
colnames(dt1)[colnames(dt1) == "P"] <- "P_moca"

colnames(dt2)[colnames(dt2) == "Estimate"] <- "est_mmse"
colnames(dt2)[colnames(dt2) == "CI_low"] <- "low_mmse"
colnames(dt2)[colnames(dt2) == "CI_up"] <- "hi_mmse"
colnames(dt2)[colnames(dt2) == "Beta(95%CI)"] <- "Beta(95%CI)_mmse"
colnames(dt2)[colnames(dt2) == "P"] <- "P_mmse"

dt1 <- cbind(dt1, dt2[, c("est_mmse", "low_mmse", "hi_mmse", "Beta(95%CI)_mmse", "P_mmse")])

fmt_beta_ci <- function(x) {
  nums <- stringr::str_extract_all(
    x, "-?\\d+\\.?\\d*(?:[eE]-?\\d+)?", simplify = TRUE
  )
  beta    <- as.numeric(nums[, 1])
  ci_low  <- as.numeric(nums[, 2])
  ci_high <- as.numeric(nums[, 3])
  
  beta_fmt <- ifelse(abs(beta) < 0.01, "<0.01", formatC(beta, format = "f", digits = 2))
  
  fix_ci <- function(v) {
    out <- formatC(v, format = "f", digits = 2)
    out[out %in% c("-0.00", "0.00")] <- "0"
    out
  }
  ci_low_fmt  <- fix_ci(ci_low)
  ci_high_fmt <- fix_ci(ci_high)
  
  pad <- function(s) ifelse(grepl("^-", s), s, paste0("\u2007", s))
  beta_fmt    <- pad(beta_fmt)
  ci_low_fmt  <- pad(ci_low_fmt)
  ci_high_fmt <- pad(ci_high_fmt)
  
  paste0(beta_fmt, " (", ci_low_fmt, ", ", ci_high_fmt, ")")
}

dt1 <- dt1 %>%
  mutate(across(c(`Beta(95%CI)_moca`, `Beta(95%CI)_mmse`), fmt_beta_ci))

add_stars <- function(p) {
  p_num <- as.numeric(p)
  ifelse(is.na(p_num), "",
         ifelse(p_num < 0.01, "**",
                ifelse(p_num < 0.05, "*", "")))
}

dt1 <- dt1 %>%
  mutate(
    star_moca = add_stars(P_moca),
    star_mmse = add_stars(P_mmse)
  ) %>%
  mutate(
    `Beta(95%CI)_moca` = paste0(`Beta(95%CI)_moca`, star_moca),
    `Beta(95%CI)_mmse` = paste0(`Beta(95%CI)_mmse`, star_mmse)
  ) %>%
  mutate(`     Beta(95%CI)` = paste0(`Beta(95%CI)_moca`, "\n", `Beta(95%CI)_mmse`))

########################
## 3) Theme
########################
alt_fill <- rep(c("#E8F0F8", "white"), length.out = nrow(dt))

font_add("Arial", regular = "arial.ttf")
showtext_auto()

tm2 <- forest_theme(
  base_size    = 10,
  base_family  = "Arial",
  refline_col  = "red",
  footnote_gp  = gpar(col = "blue"),
  legend_name  = "Cognitive score",
  legend_value = c("MoCA", "MMSE"),
  summary_col  = "#4575b4",
  core = list(
    bg_params = list(fill = alt_fill),
    txt_gp    = gpar(fontfamily = "Courier")
  ),
  ci_col = c("#4575b4", "#4daf4a"),
  ci_pch = c(16, 15),
  ci_lwd = 2,
  ci_cex = 1.5
)

########################
## 4) Plot
########################
p2 <- forest(
  dt1[, c(1, 8, 18, 10)],
  est   = list(dt1$est_moca, dt1$est_mmse),
  lower = list(dt1$low_moca, dt1$low_mmse),
  upper = list(dt1$hi_moca,  dt1$hi_mmse),
  ci_column = 4,
  ref_line  = 0,
  nudge_y   = 0.3,
  xlim      = c(-1.25, 0.5),
  ticks_at  = c(-1, -0.5, 0, 0.5),
  theme     = tm2
)

g <- edit_plot(
  p2,
  part  = "header",
  row   = 1,
  which = "background",
  gp    = gpar(fill = "#B8D8E8")
)

plot(g)

########################
## 5) Save
########################
width_inch  <- 8
height_inch <- 6

pdf(file.path(out_dir, "forest_plot.pdf"), width = width_inch, height = height_inch)
par(mar = c(5, 4, 2, 2))
plot(g)
dev.off()

png(file.path(out_dir, "forest_plot.png"),
    width = width_inch * 300, height = height_inch * 300, res = 300)
par(mar = c(5, 4, 2, 2))
plot(g)
dev.off()
