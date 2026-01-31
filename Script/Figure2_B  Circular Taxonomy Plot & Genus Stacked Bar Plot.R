###-------一.Taxa Tree plot---###########

########################
## 0) Paths & load()s
########################
rm(list=ls())

###Prepare：
# path_otu  
# path_tax 
# path_prev 
# path_enr  
# out_dir   

load(path_otu)
load(path_tax)
load(path_prev)
load(path_enr)

########################
## 1) Libraries
########################
library(MicrobiotaProcess)
library(phyloseq) 
library(ggtree)
library(dplyr)
library(reshape2)
library(ggnewscale)
library(ggplot2) 
library(stringr)
library(tibble)
library(ggtreeExtra)

########################
## 2) Data prep
########################
otu <- as.data.frame(t(otu_genus_filtered))

ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows = TRUE), 
  tax_table(as.matrix(taxa_genus_filtered))
)

df <- ps %>% as.MPSE()

taxa.tree <- df %>% 
  mp_extract_tree(type = "taxatree")

phylum_cols <- c(
  p__Actinomycetota        = "#FFFFB3",
  p__Bacillota             = "#FCCDE5",
  p__Bacteroidota          = "#BEBADA",
  p__Campylobacterota      = "#FB8072",
  p__Elusimicrobiota       = "#80B1D3",
  p__Fusobacteriota        = "#FDB462",
  p__Patescibacteria       = "#B3DE69",
  p__Pseudomonadota        = "#8DD3C7",
  p__Spirochaetota         = "#D9D9D9",
  p__Synergistota          = "#BC80BD",
  p__Thermodesulfobacteriota = "#CCEBC5"
)

########################
## 3) Plot
########################
tree_plot <- ggtree(taxa.tree, layout = "fan", open.angle = 90, linewidth = 0.5) +
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label),
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c(phylum_cols),
    labels  = str_remove(names(phylum_cols), "^p__"),
    guide = guide_legend(keywidth = 1, keyheight = 1, ncol = 2, byrow = TRUE),
    name = "Phylum"
  ) +
  new_scale_fill()

Cog <- Enrich[, "Cog", drop = FALSE]
Per <- Enrich[, "Per", drop = FALSE]

pal_cog <- c("white", "#FDD692", "#0066CC", "#E64B00")
pal_per <- c("white", "#FFE5CC", "#66B2FF", "#FF9E4D")

tree_plot2 <- 
  gheatmap(
    p              = tree_plot,
    data           = Cog,
    offset         = -0.5,
    width          = 0.2,
    colnames_angle = 0,
    color          = "grey80",
    font.size      = 2
  ) +  
  scale_fill_gradientn(
    colours = pal_cog,
    limits  = c(1, 4),
    breaks  = 2:4,
    labels  = c("Normal Cognition",
                "Mild Cognitive Impairment",
                "Dementia"),
    name    = "Enrich Group",
    guide   = "legend"
  ) + 
  new_scale_fill()

tree_plot3 <- 
  gheatmap(
    p              = tree_plot2,
    data           = Per,
    offset         = 0.9,
    width          = 0.2,
    colnames_angle = 0,
    color          = "grey80",
    font.size      = 2
  ) +  
  scale_fill_gradientn(
    colours = pal_per,
    limits  = c(1, 4),
    breaks  = 2:4,
    labels  = c("Mild Periodontitis",
                "Moderate Periodontitis",
                "Severe Periodontitis"),
    name    = "Enrich Group",
    guide   = "legend"
  )

prev_df <- Prevalence %>%                 
  as.data.frame() %>%                     
  rownames_to_column("ASV") %>%           
  rename(Prevalence = 2)

tree_plot3 <- tree_plot3 +
  new_scale_fill() +
  geom_fruit(
    data        = prev_df,
    geom        = geom_bar,
    mapping     = aes(
      y = ASV,
      x = Prevalence
    ),
    orientation = "y",
    stat        = "identity",
    pwidth      = 0.3,
    offset      = 0.49,
    fill        = "#dadbdb",
    colour      = NA
  )

print(tree_plot3)

########################
## 4) Save
########################
pdf(file = file.path(out_dir, "tree_plot3.pdf"),
    width = 8, height = 8, useDingbats = FALSE)
print(tree_plot3)
dev.off()

ggplot2::ggsave(
  filename = file.path(out_dir, "tree_plot3.png"),
  plot     = tree_plot3,
  width    = 8,
  height   = 8,
  units    = "in",
  dpi      = 500
)

#####################--------- 二. Genus bar plot--------------#####################

########################
## 0) Paths & load()s
########################
rm(list=ls())

###Prepare：
# path_otu  
# path_tax  
# path_meta 
# out_dir   

load(path_otu)
load(path_tax)
load(path_meta)

########################
## 1) Libraries
########################
library(MicrobiotaProcess)
library(phyloseq) 
library(ggtree)
library(dplyr)
library(reshape2)
library(ggnewscale)
library(ggplot2)

########################
## 2) Data prep
########################
otu <- as.data.frame(t(otu_genus_filtered))

sample <- metadata[!is.na(metadata$Diagnosis), ]

OTU <- as.data.frame(t(otu))
OTU$group <- sample$Diagnosis[match(rownames(OTU), rownames(sample))]
OTU <- OTU[, c(109, 1:108)]

data_sum <- OTU %>% 
  group_by(group) %>% 
  summarise_at(vars(1:108), sum) %>% 
  as.data.frame()

rownames(data_sum) <- data_sum$group
RelAbundance <- apply(data_sum[-1], 1, function(x) x / sum(x))
data_RelAbundance <- melt(RelAbundance, varnames = c("OTU", "group"))

data_sum2 <- as.data.frame(t(data_sum[-1]))
data_sum2$OTU <- rownames(data_sum2)

tax <- taxa_genus_filtered
tax$OTU <- rownames(tax)

data_tax_OTUsum <- merge(tax, data_sum2, by = "OTU")

data_tax_OTUsum2 <- data_tax_OTUsum[c(7, 8:10)] %>% 
  group_by(Genus) %>% 
  summarise_at(vars(1:3), sum) %>% 
  as.data.frame()

rownames(data_tax_OTUsum2) <- data_tax_OTUsum2$Genus
RelAbundance_genus <- apply(data_tax_OTUsum2[-1], 2, function(x) x / sum(x))

data_RelAbundance_genus <- melt(RelAbundance_genus, varnames = c("Genus", "group"))
data_RelAbundance_genus <- data_RelAbundance_genus %>% 
  mutate(group = factor(group,
                        levels = c("N", "MCI", "D"),
                        labels = c("CN", "MCI", "D")))

########################
## 3) Data processing: Top Genus and Prevalence
########################
top9_other <- data_RelAbundance_genus %>% 
  group_by(group) %>% 
  mutate(rank = dense_rank(desc(value)),
         Genus = ifelse(rank <= 9, as.character(Genus), "Other")) %>% 
  group_by(group, Genus) %>% 
  summarise(value = sum(value), .groups = "drop") %>% 
  group_by(group) %>% 
  arrange(desc(value), .by_group = TRUE) %>% 
  ungroup()

level_order <- top9_other %>% 
  group_by(Genus) %>% 
  summarise(total = sum(value), .groups = "drop") %>% 
  arrange(total) %>% 
  pull(Genus) %>% 
  { c("Other", setdiff(., "Other")) }

top9_other <- top9_other %>% 
  mutate(Genus = factor(Genus, levels = level_order))

cols <- c(
  Alloprevotella = "#3C5488FF",
  Fusobacterium  = "#4DBBD5FF",
  Haemophilus    = "#FDB462",
  Neisseria      = "#E64B35FF",
  Other          = "#E0E3DA",
  Porphyromonas  = "#f199bc",
  Prevotella     = "#00A087FF",
  Rothia         = "#B3DE69",
  Streptococcus  = "#BC80BD",
  Veillonella    = "#91D1C2FF"
)

########################
## 4) Plotting: Genus Bar Plot
########################
Genus_bar <- top9_other %>% 
  ggplot(aes(x = group, y = 100 * value, fill = Genus)) +
  geom_col(position = "stack", width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative Abundance (%)", fill = "Genus") +
  guides(fill = guide_legend(
    keywidth  = 1.5,
    keyheight = 1.5,
    ncol      = 2,
    byrow     = TRUE
  )) +
  scale_fill_manual(
    values = cols,
    breaks = sort(names(cols))
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.01),
    axis.line    = element_line(colour = "black", size = 0.01),
    panel.grid   = element_blank(),
    axis.text    = element_text(size = 9),
    axis.title   = element_text(size = 12)
  )

print(Genus_bar)

########################
## 5) Save Results
########################
pdf(file = file.path(out_dir, "Genus_bar.pdf"),
    width = 8, height = 6, useDingbats = FALSE)
print(Genus_bar)
dev.off()

ggplot2::ggsave(
  filename = file.path(out_dir, "Genus_bar.png"),
  plot     = Genus_bar,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 500
)


