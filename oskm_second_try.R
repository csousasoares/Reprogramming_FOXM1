##Load Packages and Set Seed

library(Seurat)
library(sctransform)
library(tidyr)
library(dplyr)
library(patchwork)
library(sctransform)
library(ggplot2)
library(UCell)
library(DESeq2)
library(MAST)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rChEA3)
library(stringr)
library(enrichplot)
library(aPEAR)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(forcats)
library(ggrepel)
library(Cairo)
library(enrichR)

set.seed(123)


## OSKM Lu et al Differential Gene Expression ----------------------------------

oskm <- readRDS("Input Files\\GSE297234_GM00731_SEVOSKM.rds")

Idents(oskm) <- "orig.ident" ## Define days as identities

DimPlot(oskm, reduction = "umap", label = T)

DefaultAssay(oskm) <- "SCT" ## Use SCT as the default assay


partial_rep_vs_fibro <- FindMarkers(
  oskm,
  slot = "data",
  ident.1 = "GM00731_D3",
  ident.2 = "GM00731_D0",
  min.pct = 0.01,
  random.seed = 123,
  logfc.threshold = 0,
  verbose = T
) ## Differential Gene Expression - Day 3 vs Day 0


partial_rep_vs_fibro$gene <- row.names(partial_rep_vs_fibro)
partial_rep_vs_fibro$p_val_adj <- p.adjust(partial_rep_vs_fibro$p_val, 
                                           method = "bonferroni")
partial_rep_vs_fibro$log10padj <- -log10(partial_rep_vs_fibro$p_val_adj)

##Define genes with 0 padjust as theoretical maximum ~ 323
partial_rep_vs_fibro$log10padj[partial_rep_vs_fibro$log10padj == Inf] <- 323

partial_rep_vs_fibro <- partial_rep_vs_fibro %>% 
  dplyr::mutate(direction = case_when(
  log10padj > 30 & avg_log2FC < -1 ~ "Downregulated",
  log10padj > 30 & avg_log2FC > 1 ~ "Upregulated",
  .default = "No Change"
)) %>% dplyr::mutate(label_oskm = case_when(
  gene %in% c("POU5F1", "SOX2", "KLF4", "MYC") ~ gene,
  .default = ""
))

direction_cols <- c(
  "Downregulated" = "blue",
  "Upregulated" = "red",
  "No Change" = "gray"
)

ggplot(partial_rep_vs_fibro, aes(x = avg_log2FC, y = log10padj, 
                                 label = label_oskm)) +
  geom_point(aes(color = direction), size = 1.2, alpha = 0.4) +
  theme_bw(base_size = 12) +
  ggrepel::geom_label_repel(box.padding = 0.8, max.overlaps = Inf, 
                            min.segment.length = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dashed") +
  scale_color_manual(values = direction_cols, name = "Direction") +
  labs(title = "Volcano Plot - Day 3 vs Day 0 Reprogramming",
       x = "Log2FC", y = "-log10(padjust)") + ylim(0, 400) + xlim(-10, 17)


ggsave(
  "day3_vs_day0_belmonte_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)


## OSKM GSEA Analysis (Biological Processes) -----------------------------------

day3_vs_day0_ranked <- partial_rep_vs_fibro %>% 
  dplyr::mutate(sign = avg_log2FC/abs(avg_log2FC)) %>% 
  dplyr::mutate(log10pval = -log10(p_val)) %>% 
  dplyr::mutate(log10pval = case_when(
    log10pval == Inf ~ 323,
    .default = log10pval
  )) %>% 
  dplyr::mutate(metric = sign*log10pval) %>% 
  dplyr::arrange(desc(metric)) ## Ranked by Sign * -log10(pval)

day3_vs_day0_ranked2 <- day3_vs_day0_ranked$metric
names(day3_vs_day0_ranked2) <- day3_vs_day0_ranked$gene
head(day3_vs_day0_ranked2, 10)

day_3_vs_0_gsea <- gseGO(geneList     = day3_vs_day0_ranked2,
                         OrgDb        = org.Hs.eg.db,
                         ont          = "BP",
                         keyType = "SYMBOL",
                         pAdjustMethod = "BH",
                         minGSSize    = 50,
                         maxGSSize    = 700,
                         eps = 1e-300,
                         pvalueCutoff = 1,
                         verbose      = TRUE,
                         seed = TRUE,
                         by = "fgsea")


df_day_3_vs_0_gsea <- as.data.frame(day_3_vs_0_gsea)





df_day_3_vs_0_gsea_plot <- df_day_3_vs_0_gsea %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_day_3_vs_0_gsea_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-2.5, 2.5) + labs(y = "Gene Sets\n", x = "\nNES",
                         title = "Partial Reprogramming - Top BPs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")

ggsave(
  "day3_vs_day0_gsea_bps.pdf",
  path = "Plots",
  width = 8,
  height = 8,
  unit = "in"
)


day_3_vs_0_gsea_005 <- day_3_vs_0_gsea %>% clusterProfiler::filter(
  p.adjust < 0.1)

day_3_vs_0_gsea_005 <- pairwise_termsim(day_3_vs_0_gsea_005)

day_3_vs_0_gsea_simple <- simplify(day_3_vs_0_gsea_005, cutoff=0.9, 
                                   by="p.adjust", select_fun=min)

df_day_3_vs_0_gsea_simple <- as.data.frame(day_3_vs_0_gsea_simple)


df_day_3_vs_0_gsea_simple_plot <- df_day_3_vs_0_gsea_simple %>% 
  dplyr::mutate(order_metric = -log10(p.adjust)*setSize) %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(order_metric, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_day_3_vs_0_gsea_simple_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-2.5, 2.5) + labs(y = "Gene Sets\n", x = "\nNES",
                         title = "Partial Reprogramming - Top BPs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


ggsave(
  "day3_vs_day0_gsea_bps_simplified.pdf",
  path = "Plots",
  width = 8,
  height = 8,
  unit = "in"
)


## 133 HDFs - Differential Gene Expression ------------------------

aging_counts <- read.csv("Input Files\\aging_counts_all.csv", row.names = 1, sep = ";")
aging_counts <- aging_counts[1:133]

metadata <- read.csv("Input Files\\sample_metadata_correct.csv", sep = ";")
##Samples like black_puerto_rican were defined as just "black"

metadata$Sample == colnames(aging_counts)

aging_counts <- aging_counts[,metadata$Sample] 

metadata$Sample == colnames(aging_counts)

metadata <- metadata %>% dplyr::mutate(age_group = case_when(
  Age <= 30 ~ "Younger",
  Age >= 70 ~ "Older",
  .default = ""
))


unique(metadata$Ethnicity)
unique(metadata$Gender)
unique(metadata$Region)

metadata_filt <- metadata %>% dplyr::filter(age_group %in% c("Younger", "Older"))
aging_counts_filt <- aging_counts[,metadata_filt$Sample]
metadata_filt$Sample == colnames(aging_counts_filt)

max(colSums(aging_counts_filt))

aging_counts_filt$gene <- row.names(aging_counts_filt)

data = aging_counts_filt[,"gene"]
data = as.vector(data)

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=data, 
                                columns="SYMBOL", keytype="ENTREZID")

aging_counts_filt$gene_symbol <- annots$SYMBOL

aging_counts_filt <- aging_counts_filt %>% na.omit()

aging_counts_filt <- aging_counts_filt[,1:84]


dds_aging <- DESeqDataSetFromMatrix(
  countData = aging_counts_filt,
  colData = metadata_filt,
  design = ~ Technology + Gender + Ethnicity + Region + age_group)

dds_aging


smallestGroupSize <- 42 ##Half of samples

keep <- rowSums(counts(dds_aging) >=10) >= smallestGroupSize
##Recomended in vignette

dds_aging <- dds_aging[keep,]

dds_aging <- DESeq(dds_aging)

res_aging <- results(dds_aging, contrast = c("age_group", "Older", "Younger"), 
                     cooksCutoff=FALSE, independentFiltering=FALSE) 

df_aging <- as.data.frame(res_aging)
df_aging$gene <- row.names(df_aging)

data = df_aging[,"gene"]
data = as.vector(data)

annots <- AnnotationDbi::select(org.Hs.eg.db, keys=data, 
                                columns="SYMBOL", keytype="ENTREZID")

df_aging$gene == annots$ENTREZID

df_aging$gene_symbol <- annots$SYMBOL


df_aging$log10padj <- -log10(df_aging$padj)
df_aging <- df_aging %>% dplyr::mutate(direction = case_when(
  log2FoldChange < -0.5 & padj < 0.05 ~ "Down",
  log2FoldChange > 0.5 & padj < 0.05 ~ "Up",
  .default = "No Change"
)) %>% dplyr::mutate(labels = case_when(
  padj < 0.05 ~ gene,
  .default = ""
))

df_aging <- df_aging %>% na.omit()


summary <- df_aging %>% group_by(direction) %>% summarise(n = n())

##Plot volcano, and change the number of up/down genes

ggplot(df_aging, aes(x = log2FoldChange, y = log10padj)) +
  ggrastr::rasterise(geom_point(aes(colour = direction), 
                                shape = 16, alpha = 0.5), dpi = 300) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c(
    "Down" = "blue",
    "No Change" = "gray",
    "Up" = "red"
  )) +
  ylim(0, 30) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  labs(title = "Volcano Plot - HDFs (>= 70y vs =< 30y)",
       y = "-log10(padjust)", x = "Log2FC") +
  annotate("text", x = -5, y = 23, 
           label = paste0("Down: ",summary[1,]$n), hjust = "left") +
  annotate("text", x = -5, y = 22, 
           label = paste0("Up: ", summary[3,]$n), hjust = "left") +
  annotate("text", x = -5, y = 21, 
           label = paste0("No Change: ", summary[2,]$n), hjust = "left") +
  annotate("rect", xmin=-5.2, xmax=-2, 
           ymin=20 , ymax=24, alpha=0.2, color="gray10", fill="white")


ggsave(
  "aging_fibro_70vs30_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)


## 133 HDFs - GSEA Biological Processes ----------------------------------------

df_aging_ranked <- df_aging %>% 
  dplyr::mutate(sign = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(log10pval = -log10(pvalue)) %>% 
  dplyr::mutate(log10pval = case_when(
    log10pval == Inf ~ 323,
    .default = log10pval
  )) %>% 
  dplyr::mutate(metric = sign*log10pval) %>% 
  dplyr::arrange(desc(metric))

df_aging_ranked2 <- df_aging_ranked$metric
names(df_aging_ranked2) <- df_aging_ranked$gene_symbol
head(df_aging_ranked2)

aging_gsea <- gseGO(geneList = df_aging_ranked2,
                    OrgDb = org.Hs.eg.db,
                    ont  = "BP",
                    keyType = "SYMBOL",
                    pAdjustMethod = "BH",
                    minGSSize    = 50,
                    maxGSSize    = 700,
                    eps = 1e-300,
                    pvalueCutoff = 1,
                    verbose      = TRUE,
                    seed = TRUE,
                    by = "fgsea",
                    nPermSimple = 20000)


df_aging_gsea <- as.data.frame(aging_gsea)



df_aging_gsea_plot <- df_aging_gsea %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_aging_gsea_plot, aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets") +
  labs(title = "Older vs Younger HDFs - Biological Processes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), name = "Direction")

ggsave(
  "aging_gsea_bps.pdf",
  path = "Plots",
  width = 8,
  height = 8,
  unit = "in"
)


aging_gsea_01 <- aging_gsea %>% clusterProfiler::filter(
  p.adjust < 0.1)

aging_gsea_01 <- pairwise_termsim(aging_gsea_01)

aging_gsea_01_simple <- simplify(aging_gsea_01, cutoff=0.9, 
                                   by="p.adjust", select_fun=min)

df_aging_gsea_01_simple <- as.data.frame(aging_gsea_01_simple)


df_aging_gsea_01_simple <- df_aging_gsea_01_simple %>% 
  dplyr::mutate(order_metric = -log10(p.adjust)*setSize) %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(order_metric, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_aging_gsea_01_simple, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets\n", x = "\nNES",
                         title = "Partial Reprogramming - Top BPs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


ggsave(
  "aging_hdfs_gsea_bps_simplified.pdf",
  path = "Plots",
  width = 8,
  height = 8,
  unit = "in"
)

df_aging_gsea_01 <- df_aging_gsea %>% 
  dplyr::filter(p.adjust < 0.1)

df_day_3_vs_0_gsea_01 <- df_day_3_vs_0_gsea %>% 
  dplyr::filter(p.adjust < 0.1)

common_bps <- intersect(df_aging_gsea_01$ID, df_day_3_vs_0_gsea_01$ID)

df_aging_gsea_filt <- df_aging_gsea_01 %>% dplyr::filter(ID %in% common_bps) %>% 
  dplyr::mutate(Group = "Aging")

df_day_3_vs_0_filt <- df_day_3_vs_0_gsea_01 %>% dplyr::filter(ID %in% common_bps) %>% 
  dplyr::mutate(Group = "Reprogramming")

aging_rep_gsea <- rbind(df_aging_gsea_filt, df_day_3_vs_0_filt)

aging_rep_gsea_2 <- aging_rep_gsea %>% dplyr::select(Description, NES, Group, p.adjust)

aging_rep_gsea_2 <- aging_rep_gsea_2 %>% 
  pivot_wider(names_from = Group, values_from = c(NES, p.adjust)) %>% 
  dplyr::mutate(Label = case_when(
    Description == "DNA repair" ~ Description,
    Description == "chromosome organization" ~ Description,
    Description == "canonical Wnt signaling pathway" ~ Description,
    Description == "mesenchyme development" ~ Description,
    Description == "telomere maintenance" ~ Description,
    Description == "rRNA processing" ~ Description,
    .default = ""
  ))



ggplot(aging_rep_gsea_2, 
       aes(x = NES_Aging, y = NES_Reprogramming, label = Label)) + 
  geom_text(aes(label = "\u2BCA", colour = NES_Aging, alpha = p.adjust_Aging), 
            size = 6, family = "Segoe UI Symbol",
            key_glyph = "rect") +
  geom_text(aes(label = "\u2BCB", 
                colour = NES_Reprogramming, alpha = p.adjust_Reprogramming), 
            size = 6, family = "Segoe UI Symbol") +
  geom_label_repel(max.overlaps = Inf, min.segment.length = unit(0, 'lines'),
                   force = 10,
                   force_pull = 10, nudge_x = 0.8,
                   nudge_y = 0.8) +
  theme_bw(base_size = 16) +
  scale_size_continuous(range = c(6, 1)) +
  scale_alpha_continuous(range = c(1, 0.2), transform = "log10") +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) 


ggsave(
  "correlation_gsea_bp_aging_reprog_2_2.pdf",
  path = "Plots",
  height = 7,
  width = 9,
  unit = "in",
  device = cairo_pdf
)

## Network of BP in Partial Reprogramming Altered in Aging --------------------


day_3_vs_0_gsea_01_in_aging <- day_3_vs_0_gsea_01 %>% 
  clusterProfiler::filter(ID %in% common_bps)

df_day_3_vs_0_gsea_005_in_aging <- as.data.frame(day_3_vs_0_gsea_01_in_aging)

network_data <- enrichmentNetwork(day_3_vs_0_gsea_01_in_aging@result, 
                                  verbose = T, 
                                  drawEllipses = TRUE, 
                                  fontSize = 2.5, 
                                  colorType = "nes",
                                  clustNameMethod = "pagerank",
                                  minClusterSize = 6,
                                  outerCutoff = 0.3) 

network_data + scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                                     midpoint = 0)

ggsave(
  "gsea_reprogramming_aging_network_bp_alt.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)
