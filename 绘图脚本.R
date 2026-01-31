#####################--------- 一. Tree 图--------------#####################
rm(list=ls())
library(MicrobiotaProcess)
library(phyloseq) 
library(ggtree)
library(dplyr)
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package
library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'
library(ggplot2) 
library(stringr)   # 用于正则替换，若已加载可省略
library(tibble)
library(ggtreeExtra)
load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/otu_genus_filtered.RData")
otu <- as.data.frame(t(otu_genus_filtered))
load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/taxa_genus_filtered.RData")
####### 一.内层树
###构建tree
##这里借助MicrobiotaProcess包进行构建
# 利用phyloseq包重新构造可转换为分析的数据格式
ps <- phyloseq(
  otu_table(as.matrix(otu), taxa_are_rows=TRUE), 
  tax_table(as.matrix(taxa_genus_filtered)))
#转换数据格式
df <- ps %>% as.MPSE()
df
# tree数据提取
taxa.tree <- df %>% 
  mp_extract_tree(type="taxatree")
taxa.tree
# 自定义离散调色盘：Phylum → 颜色
phylum_cols <- c(
  p__Actinomycetota        = "#FFFFB3",  # 青绿色
  p__Bacillota             = "#FCCDE5",  # 浅黄
  p__Bacteroidota          = "#BEBADA",  # 淡紫
  p__Campylobacterota      = "#FB8072",  # 珊瑚橘
  p__Elusimicrobiota       = "#80B1D3",  # 天蓝
  p__Fusobacteriota        = "#FDB462",  # 蜜橙
  p__Patescibacteria       = "#B3DE69",  # 草绿
  p__Pseudomonadota        = "#8DD3C7",  # 粉红
  p__Spirochaetota         = "#D9D9D9",  # 中性灰
  p__Synergistota          = "#BC80BD",  # 洋红
  p__Thermodesulfobacteriota = "#CCEBC5" # 淡薄荷绿
)


# ##二.外层热图圈层
# #1.流行率
# library(readxl)
# rm(list=ls())
# load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/otu_genus_filtered.RData")
# EukRA_wide <- as.data.frame(t(otu_genus_filtered))
# load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/taxa_genus_filtered.RData")
# all.equal(rownames(EukRA_wide),rownames(taxa_genus_filtered))
# ## 将 taxa_species_filtered 的第 2 至 7 列添加到 EukRA_wide 中
# EukRA_wide <- cbind(EukRA_wide, taxa_genus_filtered[, 2:6])
# ## 将这 6 个新列移到前面
# EukRA_wide <- EukRA_wide %>%
#   dplyr::select(colnames(taxa_genus_filtered)[2:6], everything())
# ##EukRA_wide$Genus <- paste(EukRA_wide$Family, EukRA_wide$Genus, sep = " ")    因为是Genus不太会有重复
# colnames(EukRA_wide)[1:5] <- c("phylum", "class", "order", "family", "genus")
# get_prevalence <- function(x){
#   library(tidyverse)
#   rownames(x) <- x[,1]
#   x <- x[, -1]
#   binary_matrix <- ifelse(x > 0, 1, 0)
#   prevalence <- rowSums(binary_matrix) / ncol(binary_matrix)
#   result_df <- data.frame(
#     Variable = rownames(x),
#     Prevalence = prevalence
#   ) %>%  mutate(inclusion = ifelse(Prevalence > 0.1 ,1,0))
#   return(result_df)
# }
# EukRA_wide[,6:1162] <- apply(EukRA_wide[,6:1162],2,as.numeric)   #变为数值型
# Euk_prevalence <- get_prevalence(EukRA_wide[,-c(1:4)])            #计算流行率
# ## 把行名(ASV)搬到列，再只保留 ASV 和 genus
# asv_key <- EukRA_wide %>%
#   rownames_to_column("ASV") %>%      # 把 rownames 变成列
#   select(ASV, genus)                 # 只留这两列
# ## 按属名对齐并合并到 Euk_prevalence
# Euk_prevalence <- Euk_prevalence %>%
#   left_join(asv_key,                 # 带 ASV 的键表
#             by = c(Variable = "genus"))
# rownames(Euk_prevalence) <- Euk_prevalence$ASV
# Prevalence <- as.matrix(Euk_prevalence[, 2, drop = FALSE])
# save(Prevalence,file = "C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/Prevalence.RData")
load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/Prevalence.RData")

# ##2.富集情况
# rm(list=ls())
# load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/taxa_genus_filtered.RData")
# 认知lefse <- read_excel("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/二.菌群丰度/lefse组间差异分析/认知lefse.xlsx",sheet = 2)
# 牙周炎lefse <- read_excel("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/二.菌群丰度/lefse组间差异分析/牙周炎lefse.xlsx",sheet = 2)
# cog_lefse <- 认知lefse %>%
#   select(Genus, Group, LDA) %>%            # 选择指定列
#   rename(
#     cog_group = Group,
#     cog_LDA = LDA  )
# Per_lefse <- 牙周炎lefse %>%
#   select(Genus, Group, LDA) %>%            # 选择指定列
#   rename(
#     Per_group = Group,
#     Per_LDA = LDA)
# combined_lefse <- full_join(               #全合并lefse
#   cog_lefse,
#   Per_lefse,
#   by = "Genus")
# Lefse <-combined_lefse[,c(1,2,4)]
# taxa_genus_filtered$ASV <- rownames(taxa_genus_filtered)
# taxa <- taxa_genus_filtered[,c(7,6)]
# #拼接过去
# taxa <- taxa %>%                      # 原 108×2 数据框
#   left_join(
#     Lefse %>%                         # 71×3
#       select(Genus, cog_group, Per_group),  # 只保留要并过去的列
#     by = "Genus"
#   ) %>%
#   rename(                             # 改成你想要的列名
#     Cog = cog_group,
#     Per = Per_group
#   )
# Enrich <- taxa
# #设置分类水平
# Enrich <- Enrich %>%
#   mutate(
#     Cog = case_when(
#       is.na(Cog)     ~ 1L,
#       Cog == "N"     ~ 2L,
#       Cog == "MCI"   ~ 3L,
#       Cog == "D"     ~ 4L,
#       TRUE           ~ NA_integer_  # 若出现其他未知值
#     ),
#     Per = case_when(
#       is.na(Per)         ~ 1L,
#       Per == "Mild"      ~ 2L,
#       Per == "Moderate"  ~ 3L,
#       Per == "Severe"    ~ 4L,
#       TRUE               ~ NA_integer_
#     )
#   );
# rownames(Enrich) <- Enrich$ASV
# Enrich <- as.matrix(Enrich[, c("Cog","Per"), drop = FALSE]);save(Enrich, file = "C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/Enrich.RData")




load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/Enrich.RData")



# 内层进化树图
# 基础树图构建（带 tip label 和高亮 Phylum）
tree_plot <- ggtree(taxa.tree, layout = "fan", open.angle = 90, linewidth = 0.5) +
  #geom_tiplab(size = 2.5, align = FALSE, offset = 5.5) +    （不显示最外层ASV）
  geom_hilight(
    data = td_filter(nodeClass == "Phylum"),
    mapping = aes(node = node, fill = label),
    alpha = 0.5
  ) +
  scale_fill_manual(
    values = c(phylum_cols),
    labels  = str_remove(names(phylum_cols), "^p__"),  # ← 去掉 “p__” 前缀
    guide = guide_legend(keywidth = 1, keyheight = 1,ncol = 2, byrow = TRUE),
    name = "Phylum"
  ) +
  new_scale_fill()  # 重置颜色，为热图准备新填色通道



# # 添加环形热图Prevalence
# tree_plot2 <- gheatmap(
#   p = tree_plot,
#   data = Prevalence,         # 行名应与 taxa.tree$tip.label 对应
#   offset = -0.1,               # 控制热图与树的距离
#   width = 0.1,                 # 热图宽度
#   colnames_angle = 0,          # 热图列名角度
#   high = "#E71D36",            # 高值色
#   low = "#FFF7F3",             # 低值色
#   color = "grey80",            # 分割线颜色
#   legend_title = "Prevalence(%)",  # 热图图例标题
#   font.size = 2              # 样本名字体大小
# )


# 添加环形热图Enrichi
Cog <- Enrich[, "Cog", drop = FALSE]   # 第一圈（认知）
Per <- Enrich[, "Per", drop = FALSE]   # 第二圈（牙周）
pal_cog <- c("white", "#FDD692", "#0066CC", "#E64B00")   # 橙→红  (认知)
pal_per <- c("white", "#FFE5CC", "#66B2FF","#FF9E4D" )   # 浅蓝→深蓝 (牙周)
#第一圈Cog
tree_plot2 <- 
  gheatmap(
    p              = tree_plot,
    data           = Cog,
    offset         = -0.5,   # 与树距离
    width          = 0.2,    # 热图宽度
    colnames_angle = 0,
    color          = "grey80",
    font.size      = 2
  ) +  scale_fill_gradientn(
    colours = pal_cog,
    limits  = c(1, 4),
    breaks  = 2:4,
    labels  = c("Normal Cognition",
                "Mild Cognitive Impairment",
                "Dementia"),
    name    = "Enrich Group",
    guide   = "legend"
  ) + new_scale_fill()

#第二圈Per
tree_plot3<- 
  gheatmap(
    p              = tree_plot2,
    data           = Per,
    offset         = 0.9,   # 与树距离
    width          = 0.2,    # 热图宽度
    colnames_angle = 0,
    color          = "grey80",
    font.size      = 2
  ) +  scale_fill_gradientn(
    colours = pal_per,
    limits  = c(1, 4),
    breaks  = 2:4,
    labels  = c("Mild Periodontitis",
                "Moderate Periodontitis",
                "Severe Periodontitis"),
    name    = "Enrich Group",
    guide   = "legend"
  ) 
#添加流行率柱状
## 1. 预处理流行率矩阵 —— 仅两列 ASV 与 Prevalence
prev_df <- Prevalence %>%                 # 108×1 矩阵
  as.data.frame() %>%                     # 转 data.frame
  rownames_to_column("ASV") %>%           # 行名搬到列
  rename(Prevalence = 2)
## 2. 在已有 tree_plot2 上追加灰色柱状环
tree_plot3 <- tree_plot3 +                # tree_plot2 是你已有的 Enrich 热图
  new_scale_fill() +                      # 开新通道，避免与前面渐变冲突
  geom_fruit(
    data        = prev_df,
    geom        = geom_bar,
    mapping     = aes(
      y = ASV,              # 行名必须与树 tip.label 一致
      x = Prevalence        # 柱长 = 流行率
    ),
    orientation = "y",
    stat        = "identity",
    pwidth      = 0.3,      # 环宽
    offset      = 0.49,     # 与上一圈距离
    fill        = "#dadbdb", # **固定灰色，不再映射 fill**
    colour      = NA        # 去掉柱边框，如需边框可改 "grey50"
  )
print(tree_plot3)

path <- "C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/"
## —— 保存 PDF（矢量）——
pdf(file = file.path(path, "tree_plot3.pdf"),
    width = 8, height = 8, useDingbats = FALSE)  # 8×8 英寸，可根据需要调整
print(tree_plot3)
dev.off()
## —— 保存 PNG（位图）——
ggplot2::ggsave(
  filename = file.path(path, "tree_plot3.png"),
  plot     = tree_plot3,
  width    = 8,          # 同样 8×8 英寸
  height   = 8,
  units    = "in",
  dpi      = 500         # 分辨率，杂志常用 300dpi
)






#####################--------- 二. 堆叠图--------------#####################



rm(list=ls())
library(MicrobiotaProcess)
library(phyloseq) 
library(ggtree)
library(dplyr)
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package
library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'
library(ggplot2) 

load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/otu_genus_filtered.RData")
otu <- as.data.frame(t(otu_genus_filtered))
load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/taxa_genus_filtered.RData")
load("C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/二.微生物部分/filter后/metadata(和otu表保持一致).RData")
sample<-metadata[!is.na(metadata$Diagnosis),]                                   ###去除metadata缺失值


#otu 行名为样本  组名+列为OTU
OTU <- as.data.frame(t(otu))
OTU$group <- sample$Diagnosis[ match(rownames(OTU), rownames(sample)) ];OTU <- OTU[,c(109,1:108)]
#按照sample中的group分组信息求和(也可以计算均值)
data_sum <- OTU%>% 
  group_by(group) %>% 
  summarise_at(vars(1:108), sum) %>% 
  as.data.frame()
#计算相对丰度
rownames(data_sum) <- data_sum$group
RelAbundance <- apply(data_sum[-1],1,function(x) x/sum(x))
#转换为长数据格式
data_RelAbundance <- melt(RelAbundance, varnames = c("OTU","group"))


## 计算各组样本中的物种组成情况（class水平为例）
# 将物种分类信息与OTU信息合并
data_sum2 <- as.data.frame(t(data_sum[-1]))
data_sum2$OTU <- rownames(data_sum2)
tax <- taxa_genus_filtered
tax$OTU <- rownames(tax)
data_tax_OTUsum <- merge(tax, data_sum2, by = "OTU")
# 合并重复项
data_tax_OTUsum2 <- data_tax_OTUsum[c(7,8:10)] %>% 
  group_by(Genus) %>% 
  summarise_at(vars(1:3), sum) %>% 
  as.data.frame()
#计算相对丰度
rownames(data_tax_OTUsum2) <- data_tax_OTUsum2$Genus
RelAbundance_genus <- apply(data_tax_OTUsum2[-1],2,function(x) x/sum(x))
#转换为长数据格式
data_RelAbundance_genus <- melt(RelAbundance_genus, varnames = c("Genus","group"))
data_RelAbundance_genus <- data_RelAbundance_genus %>% 
  mutate(group = factor(group,
                        levels = c("N",  "MCI", "D"),      # 原水平顺序
                        labels = c("CN", "MCI", "D")   # 新名字
  ))
###保留前9个相对高丰度Genus+Other
top9_other <- data_RelAbundance_genus %>% 
  group_by(group) %>% 
  mutate(rank = dense_rank(desc(value)),
         Genus = ifelse(rank <= 9, as.character(Genus), "Other")) %>% 
  group_by(group, Genus) %>% 
  summarise(value = sum(value), .groups = "drop") %>% 
  group_by(group) %>% 
  arrange(desc(value), .by_group = TRUE) %>%   # 关键行：组内按 value 降序
  ungroup()
#丰度顺序
level_order <- top9_other %>% 
  group_by(Genus) %>% 
  summarise(total = sum(value), .groups = "drop") %>% 
  arrange(total) %>%                 # 1) 丰度由小→大
  pull(Genus) %>%                    # 2) 取出 Genus 向量
  { c("Other", setdiff(., "Other"))} # 3) 把 Other 放首位
#legend顺序
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


## 绘制物种组成柱状堆积图
## 3️⃣ 画图：堆叠顺序 = level_order；图例顺序 = 字母序
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
  theme_classic() +                              # ← 去掉背景网格
  theme(
    panel.border     = element_rect(             # 面板边框
      colour = "black",
      fill   = NA,
      size   = 0.01                              # ← 细实线
    ),
    axis.line        = element_line(             # 坐标轴线
      colour = "black",
      size   = 0.01                             # ← 细实线
    ),
    panel.grid       = element_blank(),          # 确保无任何网格线
    axis.text        = element_text(size = 9),
    axis.title       = element_text(size = 12)
  )
print(Genus_bar)
path <- "C:/Users/yh/Desktop/硕士/研一冲/口腔菌群/复盘全流程/牙石和BOP重新算/图表/Tree图/"



## ——— 保存 Genus_bar 为 PDF（矢量图）———
pdf(file = file.path(path, "Genus_bar.pdf"),
    width = 8, height = 6, useDingbats = FALSE)  # 根据需要调尺寸
print(Genus_bar)
dev.off()

## ——— 保存 Genus_bar 为 PNG（位图）———
ggplot2::ggsave(
  filename = file.path(path, "Genus_bar.png"),
  plot     = Genus_bar,
  width    = 8,        # 英寸
  height   = 6,        # 若想保持与 PDF 同高可改成 8
  units    = "in",
  dpi      = 500      # 期刊常用分辨率
)

