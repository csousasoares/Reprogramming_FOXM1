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


## OSKM Lu et al 2025 Loading and DGE ------------------------------------------

oskm <- readRDS("GSE297234_GM00731_SEVOSKM.rds")

Idents(oskm) <- "orig.ident"

DefaultAssay(oskm) <- "SCT"

partial_rep_vs_fibro <- FindMarkers(
  oskm,
  slot = "data",
  ident.1 = "GM00731_D3",
  ident.2 = "GM00731_D0",
  min.pct = 0.01,
  random.seed = 123,
  logfc.threshold = 0,
  verbose = T
)


partial_rep_vs_fibro$gene <- row.names(partial_rep_vs_fibro)
partial_rep_vs_fibro$p_val_adj <- p.adjust(partial_rep_vs_fibro$p_val,
                                           method = "bonferroni")
partial_rep_vs_fibro$log10padj <- -log10(partial_rep_vs_fibro$p_val_adj)


partial_rep_vs_fibro$log10padj[partial_rep_vs_fibro$log10padj == Inf] <- 323

partial_rep_vs_fibro <- partial_rep_vs_fibro %>% 
  dplyr::mutate(direction = case_when(
  log10padj > 30 & avg_log2FC < -1 ~ "Down",
  log10padj > 30 & avg_log2FC > 1 ~ "Up",
  .default = "No Change"
)) %>% dplyr::mutate(label_oskm = case_when(
  gene %in% c("POU5F1", "SOX2", "KLF4", "MYC") ~ gene,
  .default = ""
))

direction_cols <- c(
  "Down" = "blue",
  "Up" = "red",
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


## 133 HDFs Data Loading and DGE -----------------------------------------------


aging_counts <- read.csv("aging_counts_all.csv", row.names = 1, sep = ";")
aging_counts <- aging_counts[1:133]

metadata <- read.csv("sample_metadata_correct.csv", sep = ";")
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

metadata_filt <- metadata %>% 
  dplyr::filter(age_group %in% c("Younger", "Older"))
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

res_aging <- results(dds_aging, 
                     contrast = c("age_group", "Older", "Younger"), 
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
  annotate("text", x = -5, y = 23, label = paste0("Down: ",summary[1,]$n), 
           hjust = "left") +
  annotate("text", x = -5, y = 22, label = paste0("Up: ", summary[3,]$n), 
           hjust = "left") +
  annotate("text", x = -5, y = 21, label = paste0("No Change: ", 
           summary[2,]$n), 
           hjust = "left") +
  annotate("rect", xmin=-5.2, xmax=-2, ymin=20 , ymax=24, alpha=0.2, 
           color="gray10", fill="white")


ggsave(
  "aging_fibro_70vs30_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)


## GSEA OSKM scRNA-seq ---------------------------------------------------------

day3_vs_day0_ranked <- partial_rep_vs_fibro %>% 
  dplyr::mutate(sign = avg_log2FC/abs(avg_log2FC)) %>% 
  dplyr::mutate(log10pval = -log10(p_val)) %>% 
  dplyr::mutate(log10pval = case_when(
    log10pval == Inf ~ 323,
    .default = log10pval
  )) %>% 
  dplyr::mutate(metric = sign*log10pval) %>% 
  dplyr::arrange(desc(metric))

day3_vs_day0_ranked2 <- day3_vs_day0_ranked$metric
names(day3_vs_day0_ranked2) <- day3_vs_day0_ranked$gene
head(day3_vs_day0_ranked2)

day_3_vs_0_gsea <- gseGO(geneList     = day3_vs_day0_ranked2,
                         OrgDb        = org.Hs.eg.db,
                         ont          = "BP",
                         keyType = "SYMBOL",
                         pAdjustMethod = "BH",
                         minGSSize    = 30,
                         maxGSSize    = 700,
                         eps = 1e-300,
                         pvalueCutoff = 1,
                         verbose      = TRUE,
                         seed = TRUE,
                         by = "fgsea")


df_day_3_vs_0_gsea <- as.data.frame(day_3_vs_0_gsea)

nsig <- df_day_3_vs_0_gsea %>% dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::summarise(n = n()) %>% dplyr::pull(n)

day_3_vs_0_gsea <- pairwise_termsim(day_3_vs_0_gsea)
day_3_vs_0_gsea_005 <- day_3_vs_0_gsea %>% 
  clusterProfiler::filter(p.adjust < 0.05)


df_day_3_vs_0_gsea_plot <- df_day_3_vs_0_gsea %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 40) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_day_3_vs_0_gsea_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-2.5, 2.5) + labs(y = "Gene Sets") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")



df_day_3_vs_0_gsea_005 <- df_day_3_vs_0_gsea %>% dplyr::filter(p.adjust < 0.05)


## Aging GSEA and Common BPs----------------------------------------------------


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
                    minGSSize    = 30,
                    maxGSSize    = 700,
                    eps = 1e-300,
                    pvalueCutoff = 1,
                    verbose      = TRUE,
                    seed = TRUE,
                    by = "fgsea",
                    nPermSimple = 20000)


df_aging_gsea <- as.data.frame(aging_gsea)

nsig <- df_aging_gsea %>% dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::summarise(n = n()) %>% dplyr::pull(n)

aging_gsea <- pairwise_termsim(aging_gsea)
aging_gsea_005 <- aging_gsea %>% clusterProfiler::filter(p.adjust < 0.05)


df_aging_gsea_plot <- df_aging_gsea %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_aging_gsea_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")



df_aging_gsea_005 <- df_aging_gsea %>% dplyr::filter(p.adjust < 0.05)


common_bps <- intersect(df_aging_gsea_005$ID, df_day_3_vs_0_gsea_005$ID)

df_aging_gsea_005_filt <- df_aging_gsea_005 %>% 
  dplyr::filter(ID %in% common_bps) %>% 
  dplyr::mutate(Group = "Aging")

df_day_3_vs_0_005_filt <- df_day_3_vs_0_gsea_005 %>% 
  dplyr::filter(ID %in% common_bps) %>% 
  dplyr::mutate(Group = "Reprogramming")

aging_rep_gsea <- rbind(df_day_3_vs_0_005_filt, df_aging_gsea_005_filt)

aging_rep_gsea_2 <- aging_rep_gsea %>% 
  dplyr::select(Description, NES, Group, p.adjust)

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
  geom_text(aes(label = "\u2BCB", colour = NES_Reprogramming, 
                alpha = p.adjust_Reprogramming), 
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
  "correlation_gsea_bp_aging_reprog_2.pdf",
  path = "Plots",
  height = 7,
  width = 9,
  unit = "in",
  device = cairo_pdf
)

## Reprogramming + Aging Network and TF Enrichment -----------------------------

day_3_vs_0_gsea_005_in_aging <- day_3_vs_0_gsea_005 %>% 
  clusterProfiler::filter(ID %in% common_bps)

df_day_3_vs_0_gsea_005_in_aging <- as.data.frame(day_3_vs_0_gsea_005_in_aging)

network_data <- enrichmentNetwork(day_3_vs_0_gsea_005_in_aging@result, 
                                  verbose = T, 
                                  drawEllipses = TRUE, 
                                  fontSize = 2.5, 
                                  colorType = "nes",
                                  clustNameMethod = "pagerank",
                                  minClusterSize = 5,
                                  outerCutoff = 0.2) 

network_data + scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                                     midpoint = 0)

ggsave(
  "gsea_reprogramming_aging_network_bp.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)

clusters_of_bps <- network_data$plot_env$clust
clusters_of_bps <- as.data.frame(clusters_of_bps)
clusters_of_bps$Description <- row.names(clusters_of_bps)

clusters_of_bps <- clusters_of_bps %>% 
  dplyr::left_join(df_day_3_vs_0_gsea_005_in_aging, join_by(Description))

colnames(clusters_of_bps)
new_names <- clusters_of_bps %>% dplyr::group_by(clusters_of_bps) %>% 
  dplyr::slice_max(setSize, n = 1)
new_names

rm(oskm)


clusters_plot <- clusters_of_bps %>% 
  dplyr::group_by(clusters_of_bps) %>% dplyr::slice_min(p.adjust, n = 5)


clusters_plot



my_list <- list()

bp_list <- clusters_plot$Description

for (bp in bp_list) {
  
  set <- clusters_plot %>% dplyr::filter(Description == bp)
  genes <- str_split(set$core_enrichment, "/")
  genes <- genes[[1]]
  print(set$Description)
  print(genes)
  results <- queryChEA3(genes)
  results2 <- results[[1]]
  results2$n_genes <- str_count(results2$Overlapping_Genes, pattern = ",") + 1
  TFs <- results2 %>% dplyr::arrange(desc(n_genes)) %>% 
    dplyr::slice_max(n_genes, n = 5)
  TFs_2 <-  paste0(TFs$TF[1:5], collapse=",")
  TFs_3 <- data.frame("TFs" = TFs_2, "Description" = bp)
  my_list[[bp]] <- TFs_3
}



tf_results <- do.call(rbind.data.frame, my_list)

clusters_plot <- clusters_plot %>% 
  dplyr::left_join(tf_results, join_by(Description))


clusters_plot$Description <- str_wrap(clusters_plot$Description, width = 20)



colnames(clusters_plot)
clusters_plot_positive <- clusters_plot %>% dplyr::filter(NES > 0)

ggplot(clusters_plot_positive, 
       aes(x = setSize, y = fct_reorder(Description, setSize))) + 
  geom_col(aes(fill = p.adjust), color = "gray", width = 0.6) + 
  facet_wrap(.~clusters_of_bps, scales = "free") +
  geom_text(aes(label = TFs), 
            x = 0,
            hjust = -0.05,
            size = 3.5) + 
  coord_cartesian(clip = "on") +
  theme_bw(base_size = 14) +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5)) +
  scale_fill_gradient(high = "white", low = "orange", 
                      transform = "log10", name = "padjust") + 
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black")
  ) +
  labs(x = "Gene Count", y = "Gene Sets",
       title = "TF Enrichment in Partial Reprogramming and Aging",
       subtitle = "GSEA + ChEA3")



ggsave(
  "tf_bp_barplot_reprogramming_and_aging.pdf",
  path = "Plots",
  width = 16,
  height = 11,
  unit = "in"
)


## Day 3 of Reprogramming - FOXM1 High vs Low DGE ------------------------------

Idents(oskm) <- "orig.ident"

oskm_day_3 <- subset(oskm, idents = "GM00731_D3")

DefaultAssay(oskm_day_3) <- "SCT"

VlnPlot(oskm_day_3, features = "FOXM1")

expr <- FetchData(oskm_day_3, vars = "FOXM1", slot = "data")
expr

foxm1_positive_cells <- expr %>% dplyr::filter(FOXM1 > 0) %>% row.names()
foxm1_negative_cells <- expr %>% dplyr::filter(FOXM1 == 0) %>% row.names()

oskm_day_3@meta.data[["FOXM1_Status"]] <- "Dummy"
oskm_day_3$FOXM1_Status[foxm1_positive_cells] <- "FOXM1_Pos"
oskm_day_3$FOXM1_Status[foxm1_negative_cells] <- "FOXM1_Neg"

Idents(oskm_day_3) <- "FOXM1_Status"

foxm1_pos_vs_neg <- FindMarkers(oskm_day_3,
                                slot = "data",
                                ident.1 = "FOXM1_Pos",
                                ident.2 = "FOXM1_Neg",
                                min.pct = 0.01,
                                random.seed = 123,
                                logfc.threshold = 0,
                                verbose = T)

foxm1_pos_vs_neg$gene <- row.names(foxm1_pos_vs_neg)
foxm1_pos_vs_neg$p_val_adj <- p.adjust(foxm1_pos_vs_neg$p_val, 
                                       method = "bonferroni")
foxm1_pos_vs_neg$log10padj <- -log10(foxm1_pos_vs_neg$p_val_adj)


foxm1_pos_vs_neg$log10padj[foxm1_pos_vs_neg$log10padj == Inf] <- 323

foxm1_pos_vs_neg <- foxm1_pos_vs_neg %>% dplyr::mutate(direction = case_when(
  log10padj > 3 & avg_log2FC < -0.5 ~ "Down",
  log10padj > 3 & avg_log2FC > 0.5 ~ "Up",
  .default = "No Change"
)) %>% dplyr::mutate(label_oskm = case_when(
  gene %in% c("POU5F1", "SOX2", "KLF4", "MYC") ~ gene,
  .default = ""
))

direction_cols <- c(
  "Down" = "blue",
  "Up" = "red",
  "No Change" = "gray"
)

ggplot(foxm1_pos_vs_neg, aes(x = avg_log2FC, y = log10padj)) +
  geom_point(aes(color = direction), size = 1.2, alpha = 0.4) +
  theme_bw(base_size = 12) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = 3, linetype = "dashed") +
  scale_color_manual(values = direction_cols, name = "Direction") +
  labs(title = "Volcano Plot - Day 3 vs Day 0 Reprogramming",
       x = "Log2FC", y = "-log10(padjust)") + ylim(0, 350) + xlim(-13, 13)



ggsave(
  "day3_vs_day0_belmonte_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)


## FOXM1 High vs Low GSEA BPs -------------------------------------------------- 

foxm1_pos_vs_neg <- foxm1_pos_vs_neg %>% 
  dplyr::mutate(sign = avg_log2FC/abs(avg_log2FC)) %>% 
  dplyr::mutate(log10pval = case_when(
    p_val > 0 ~ -log10(p_val),
    .default = 323
  )) %>% 
  dplyr::mutate(metric = sign * log10pval)

foxm1_ranked <- foxm1_pos_vs_neg %>% 
  dplyr::arrange(desc(metric)) 
foxm1_ranked2 <- foxm1_ranked$metric
names(foxm1_ranked2) <- foxm1_ranked$gene
head(foxm1_ranked2, 50)

foxm1_high_gsea <- gseGO(geneList = foxm1_ranked2,
                         OrgDb = org.Hs.eg.db,
                         ont  = "BP",
                         keyType = "SYMBOL",
                         pAdjustMethod = "BH",
                         minGSSize    = 30,
                         maxGSSize    = 700,
                         eps = 1e-300,
                         pvalueCutoff = 1,
                         verbose      = TRUE,
                         seed = TRUE,
                         by = "fgsea",
                         nPermSimple = 20000)


df_foxm1_high_gsea <- as.data.frame(foxm1_high_gsea)

nsig <- df_foxm1_high_gsea %>% dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::summarise(n = n()) %>% dplyr::pull(n)

foxm1_high_gsea_005 <- foxm1_high_gsea %>% 
  clusterProfiler::filter(p.adjust < 0.05)


df_foxm1_high_gsea_plot <- df_foxm1_high_gsea %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_foxm1_high_gsea_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


## FOXM1 OE vs Reprogramming ---------------------------------------------------


dndk_counts <- read.csv("Input Files\\dndk_counts.csv", sep = ";")

dndk_counts2 <- dndk_counts %>% dplyr::group_by(gene_name) %>% 
  dplyr::summarise(n = n()) %>% dplyr::filter(n == 1) %>% 
  dplyr::pull(gene_name)

dndk_counts <- dndk_counts %>% 
  dplyr::filter(gene_name %in% dndk_counts2)

row.names(dndk_counts) <- dndk_counts$gene_name

dndk_counts <- dndk_counts %>% dplyr::select(!gene_name)

dndk_counts["FOXM1",]

dndk_metadata <- read.csv("Input Files\\dndk_metadata.csv", sep = ";")


dndk_metadata$Sample == colnames(dndk_counts)


dds_dndk <- DESeqDataSetFromMatrix(
  countData = dndk_counts,
  colData = dndk_metadata,
  design = ~ Group)

dds_dndk


smallestGroupSize <- 2 ##Half of samples
keep <- rowSums(counts(dds_dndk) >=10) >= smallestGroupSize
##Recomended in vignette
dds_dndk <- dds_dndk[keep,]

dds_dndk <- DESeq(dds_dndk)

res_dndk <- results(dds_dndk, contrast = c("Group", "FOXM1_OE", "Vector"), 
                    cooksCutoff=FALSE, independentFiltering=FALSE) 

df_dndk <- as.data.frame(res_dndk)
df_dndk$gene <- row.names(res_dndk)



df_dndk$log10padj <- -log10(df_dndk$padj)
df_dndk <- df_dndk %>% dplyr::mutate(direction = case_when(
  log2FoldChange < -0.5 & padj < 0.05 ~ "Down",
  log2FoldChange > 0.5 & padj < 0.05 ~ "Up",
  .default = "No Change"
)) %>% dplyr::mutate(labels = case_when(
  padj < 0.05 ~ gene,
  .default = ""
))

df_dndk <- df_dndk %>% na.omit()


summary <- df_dndk %>% group_by(direction) %>% summarise(n = n())

##Plot volcano, and change the number of up/down genes

ggplot(df_dndk, aes(x = log2FoldChange, y = log10padj)) +
  ggrastr::rasterise(geom_point(aes(colour = direction), 
                                shape = 16, alpha = 0.5), dpi = 300) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c(
    "Down" = "blue",
    "No Change" = "gray",
    "Up" = "red"
  )) +
  ylim(0, 200) +
  xlim(-7, 7) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  labs(title = "Volcano Plot - 87y FOXM1 OE vs Vector",
       y = "-log10(padjust)", x = "Log2FC") +
  annotate("text", x = -6, y = 150, label = paste0("Down: ",summary[1,]$n), 
           hjust = "left") +
  annotate("text", x = -6, y = 140, label = paste0("Up: ", summary[3,]$n), 
           hjust = "left") +
  annotate("text", x = -6, y = 130, 
           label = paste0("No Change: ", summary[2,]$n), 
           hjust = "left") +
  annotate("rect", xmin=-6.2, xmax= -1, ymin=120 , 
           ymax=160, alpha=0.2, color="gray10", fill="white")



ggsave(
  "87y_foxm1_oe_vs_vector_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)



## GSEA BPs

dndk_ranked <- df_dndk %>% 
  dplyr::mutate(sign = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(log10pval = -log10(pvalue)) %>% 
  dplyr::mutate(metric = sign * log10pval) %>% 
  dplyr::arrange(desc(metric))

dndk_ranked2 <- dndk_ranked$metric
names(dndk_ranked2) <- dndk_ranked$gene
head(dndk_ranked2, 50)


dndk_hall <- GSEA(
  dndk_ranked2,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-30,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  gson = NULL,
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea")


df_dndk_hall <- as.data.frame(dndk_hall)


df_dndk_hall_plot <- df_dndk_hall %>% 
  dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))

df_dndk_hall_plot$Description <- gsub("HALLMARK_", "", 
                                      df_dndk_hall_plot$Description)

ggplot(df_dndk_hall_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-4, 4) + labs(y = "Hallmark Gene Sets",
                     title = "87y FOXM1 OE vs Vector",
                     subtitle = "Hallmark Collection") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


ggsave(
  "87y_foxm1_oe_vs_vector_hallmark.pdf",
  path = "Plots",
  width = 7,
  height = 8,
  unit = "in"
)



foxm1_dndk_gsea_bp <- gseGO(geneList = dndk_ranked2,
                            OrgDb = org.Hs.eg.db,
                            ont  = "BP",
                            keyType = "SYMBOL",
                            pAdjustMethod = "BH",
                            minGSSize    = 30,
                            maxGSSize    = 700,
                            eps = 1e-300,
                            pvalueCutoff = 1,
                            verbose      = TRUE,
                            seed = TRUE,
                            by = "fgsea",
                            nPermSimple = 20000)


df_foxm1_dndk_gsea_bp <- as.data.frame(foxm1_dndk_gsea_bp)

nsig <- df_foxm1_dndk_gsea_bp %>% dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::summarise(n = n()) %>% dplyr::pull(n)

foxm1_dndk_gsea_bp_005 <- foxm1_dndk_gsea_bp %>% 
  clusterProfiler::filter(p.adjust < 0.05)


df_foxm1_dndk_gsea_bp_plot <- df_foxm1_dndk_gsea_bp %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_foxm1_dndk_gsea_bp_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


foxm1_dndk_gsea_bp_005 <- pairwise_termsim(foxm1_dndk_gsea_bp_005)


network_data <- enrichmentNetwork(foxm1_dndk_gsea_bp_005@result, 
                                  verbose = T, 
                                  drawEllipses = TRUE, 
                                  fontSize = 2.5, 
                                  colorType = "nes",
                                  clustNameMethod = "pagerank",
                                  minClusterSize = 5,
                                  outerCutoff = 0.3) 

network_data + scale_color_gradient2(low = "blue", 
                                     mid = "white", 
                                     high = "red", midpoint = 0)

ggsave(
  "network_bps_87y_foxm1_oe_dndk.pdf",
  path = "Plots",
  width = 7,
  height = 7,
  unit = "in"
)



df_foxm1_dndk_gsea_bp_signifs <- df_foxm1_dndk_gsea_bp %>% 
  dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::pull(ID)

df_reprog_signif <- df_day_3_vs_0_gsea %>% 
  dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::pull(ID)

df_aging_signif <- df_aging_gsea %>% 
  dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::pull(ID)

common_reprog_dndk <- intersect(df_foxm1_dndk_gsea_bp_signifs,
                                df_reprog_signif)

common_reprog_dndk_aging <- intersect(common_reprog_dndk,
                                      df_aging_signif)



foxm1_dndk_gsea_bp_reprog_aging <- foxm1_dndk_gsea_bp %>% 
  clusterProfiler::filter(ID %in% common_reprog_dndk_aging)

network_data <- enrichmentNetwork(foxm1_dndk_gsea_bp_reprog_aging@result, 
                                  verbose = T, 
                                  drawEllipses = TRUE, 
                                  fontSize = 2.5, 
                                  colorType = "nes",
                                  clustNameMethod = "pagerank",
                                  minClusterSize = 5,
                                  outerCutoff = 0.3) 

network_data + scale_color_gradient2(low = "blue", 
                                     mid = "white", 
                                     high = "red", midpoint = 0)

ggsave(
  "network_bps_87y_foxm1_oe_affected_in_aging_reprogramming.pdf",
  path = "Plots",
  width = 7,
  height = 7,
  unit = "in"
)


## FOXM1 OE Half Circle Plot with Aging and Reprogramming -------------------------------


common_reprog_dndk <- intersect(df_foxm1_dndk_gsea_bp_signifs,
                                df_reprog_signif)


df_reprog_gsea_005_filt <- df_day_3_vs_0_gsea_005 %>% 
  dplyr::filter(ID %in% common_reprog_dndk) %>% 
  dplyr::mutate(Group = "Reprogramming")

df_dndk_005_filt <- df_foxm1_dndk_gsea_bp_005 %>% 
  dplyr::filter(ID %in% common_reprog_dndk) %>% 
  dplyr::mutate(Group = "FOXM1 OE")

dndk_rep_gsea <- rbind(df_reprog_gsea_005_filt, df_dndk_005_filt)

dndk_rep_gsea_2 <- dndk_rep_gsea %>% 
  dplyr::select(Description, NES, Group, p.adjust)

dndk_rep_gsea_2 <- dndk_rep_gsea_2 %>% 
  pivot_wider(names_from = Group, values_from = c(NES, p.adjust)) %>% 
  dplyr::mutate(Label = case_when(
    Description == "DNA repair" ~ Description,
    Description == "chromosome organization" ~ Description,
    Description == "canonical Wnt signaling pathway" ~ Description,
    Description == "mesenchyme development" ~ Description,
    Description == "telomere maintenance" ~ Description,
    Description == "rRNA processing" ~ Description,
    Description == "tube morphogenesis" ~ Description,
    Description == "hormone metabolic process" ~ Description,
    Description == "DNA-templated DNA replication" ~ Description, 
    Description == "positive regulation of cell-substrate adhesion" ~ Description,
    Description == "blood vessel morphogenesis" ~ Description,
    .default = ""
  ))

colnames(dndk_rep_gsea_2)

ggplot(dndk_rep_gsea_2, aes(x = `NES_FOXM1 OE`, 
                            y = NES_Reprogramming, label = Label)) + 
  geom_text(aes(label = "\u2BCA", colour = `NES_FOXM1 OE`, 
                alpha = `p.adjust_FOXM1 OE`), 
            size = 6, family = "Segoe UI Symbol",
            key_glyph = "rect") +
  geom_text(aes(label = "\u2BCB", colour = NES_Reprogramming, 
                alpha = p.adjust_Reprogramming), 
            size = 6, family = "Segoe UI Symbol") +
  geom_label_repel(max.overlaps = Inf, min.segment.length = unit(0, 'lines'),
                   force = 10,
                   force_pull = 10, nudge_x = 0.8,
                   nudge_y = 0.8) +
  theme_bw(base_size = 16) +
  scale_size_continuous(range = c(6, 1)) +
  scale_alpha_continuous(range = c(1, 0.2), transform = "log10") +
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        mid = "white", midpoint = 0) +
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) ## Good


ggsave(
  "correlation_gsea_bp_dndk_reprog.pdf",
  path = "Plots",
  height = 7,
  width = 9,
  unit = "in",
  device = cairo_pdf
)


## Convert to Dotplot

dndk_rep_gsea <- rbind(df_reprog_gsea_005_filt, df_dndk_005_filt)

in_both <- dndk_rep_gsea %>% dplyr::group_by(Description) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n == 2) %>% 
  dplyr::pull(Description)

dndk_rep_gsea_3 <- dndk_rep_gsea %>% 
  dplyr::select(Description, NES, Group, p.adjust) %>% 
  dplyr::filter(Description %in% in_both)
  

dndk_rep_gsea_4 <- dndk_rep_gsea_2 %>% 
  na.omit() %>% 
  dplyr::filter(Description %in% in_both) %>% 
  dplyr::mutate(direction = 
                  (NES_Reprogramming/`NES_FOXM1 OE`)/abs(
                    (NES_Reprogramming/`NES_FOXM1 OE`)  
                  )) %>% 
  dplyr::mutate(direction2 = case_when(
    direction == 1 ~ "Similar",
    direction == -1 ~ "Opposite"
  ))


dndk_rep_gsea_5 <- dndk_rep_gsea_3 %>% 
  dplyr::left_join(similar_and_opposite_go, by = "Description") %>% 
  na.omit()

dndk_rep_gsea_6_gene_sets <- dndk_rep_gsea_3 %>% 
  dplyr::left_join(similar_and_opposite_go, by = "Description") %>% 
  na.omit() %>% 
  dplyr::group_by(direction2, Group) %>% 
  dplyr::slice_min(p.adjust, n = 30) %>% 
  dplyr::ungroup() %>% 
  dplyr::pull(Description)

dndk_rep_gsea_6 <- dndk_rep_gsea_5 %>% 
  dplyr::filter(Description %in% dndk_rep_gsea_6_gene_sets)

ggplot(dndk_rep_gsea_6, 
       aes(x = Group, y = Description)) +
  geom_point(aes(fill = NES, size = p.adjust), shape = 21) +
  facet_wrap(.~direction2, scales = "free") +
  labs(title = "BPs Common to FOXM1 OE and Partial Reprogramming",
       subtitle = "In Top 30 of Each Condition") +
  scale_fill_gradient2(low = "blue",
                        mid = "white",
                        high = "red",
                        midpoint = 0) +
  scale_size_continuous(range = c(5, 2), transform = "log10") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1,
                                   vjust = 1))

ggsave(
  "foxm1_partial_reprogramming_bps_opposite_and_similar.pdf",
  path = "Plots",
  width = 12,
  height = 8,
  unit = "in"
)


## Compare Aging, FOXM1 and PR at BPs ------------------------------------------


aging_gsea <- pairwise_termsim(aging_gsea)


df_aging_gsea <- as.data.frame(aging_gsea)

sims <- aging_gsea@termsim

names <- row.names(sims)


changed_in_aging_signif <- df_aging_gsea %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::mutate(Group = "Aging") %>% 
  dplyr::filter(Description %in% names)


df_reprog_gsea_005_filt <- df_day_3_vs_0_gsea_005 %>% 
  dplyr::filter(ID %in% changed_in_aging_signif$ID) %>% 
  dplyr::mutate(Group = "Aged Cells + Reprogramming")

df_dndk_005_filt <- df_foxm1_dndk_gsea_bp_005 %>% 
  dplyr::filter(ID %in% changed_in_aging_signif$ID) %>% 
  dplyr::mutate(Group = "Aged Cells + FOXM1 OE")


reprog_or_dndk <- c(df_reprog_gsea_005_filt$ID, df_dndk_005_filt$ID) %>% 
  unique()


changed_in_aging_signif <- changed_in_aging_signif %>% 
  dplyr::filter(ID %in% reprog_or_dndk)

aging_plus_two <- rbind(changed_in_aging_signif,
                        df_reprog_gsea_005_filt,
                        df_dndk_005_filt) %>% 
  dplyr::filter(ID %in% changed_in_aging_signif$ID)


top25 <- changed_in_aging_signif %>% dplyr::slice_max(, n = 25) %>% 
  dplyr::pull(ID)
bottom25 <- changed_in_aging_signif %>% dplyr::slice_min(NES, n = 25) %>% 
  dplyr::pull(ID)

top50 <- c(top25, bottom25)

top50 <- changed_in_aging_signif %>% dplyr::slice_max(setSize, n = 30) %>% 
  dplyr::pull(ID)

aging_plus_two_2 <- aging_plus_two %>% dplyr::filter(ID %in% top50) %>% 
  dplyr::mutate(NES_dummy = case_when(
    Group == "Aging" ~ NES,
    .default = 0
  ))

ggplot(aging_plus_two_2, 
       aes(x = Group, y = forcats::fct_reorder(Description, -NES_dummy))) +
  geom_point(aes(fill = NES, size = setSize), shape = 21) +
  scale_fill_gradient2(low = "blue", 
                       mid = "white", 
                       high = "red", midpoint = 0) +
  theme_bw(base_size = 14) +
  labs(title = "Aging, Reprogramming and FOXM1 OE",
       x = "\nComparison", y = "Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  "aging_foxm1_oe_pr_bps_plot.pdf",
  path = "Plots",
  width = 11,
  height = 8,
  unit = "in"
)


## siFOXM1 Differential Gene Expression ----------------------------------------


sifoxm1_counts <- read.csv("Input Files\\sifoxm1_counts.csv", sep = ";")

sifoxm1_counts2 <- sifoxm1_counts %>% dplyr::group_by(gene_name) %>% 
  dplyr::summarise(n = n()) %>% dplyr::filter(n == 1) %>% 
  dplyr::pull(gene_name)

sifoxm1_counts <- sifoxm1_counts %>% 
  dplyr::filter(gene_name %in% sifoxm1_counts2)

row.names(sifoxm1_counts) <- sifoxm1_counts$gene_name

sifoxm1_counts <- sifoxm1_counts %>% dplyr::select(!gene_name)

sifoxm1_counts["FOXM1",]

sifoxm1_metadata <- read.csv("Input Files\\sifoxm1_metadata.csv", sep = ";")


sifoxm1_metadata$Sample == colnames(sifoxm1_counts)


dds_sifoxm1 <- DESeqDataSetFromMatrix(
  countData = sifoxm1_counts,
  colData = sifoxm1_metadata,
  design = ~ Group)

dds_sifoxm1


smallestGroupSize <- 4 ##Half of samples
keep <- rowSums(counts(dds_sifoxm1) >=10) >= smallestGroupSize
##Recomended in vignette
dds_sifoxm1 <- dds_sifoxm1[keep,]

dds_sifoxm1 <- DESeq(dds_sifoxm1)

res_sifoxm1 <- results(dds_sifoxm1, contrast = c("Group", "siFOXM1", "Mock"), 
                    cooksCutoff=FALSE, independentFiltering=FALSE) 

df_sifoxm1 <- as.data.frame(res_sifoxm1)
df_sifoxm1$gene <- row.names(res_sifoxm1)



df_sifoxm1$log10padj <- -log10(df_sifoxm1$padj)
df_sifoxm1 <- df_sifoxm1 %>% dplyr::mutate(direction = case_when(
  log2FoldChange < -0.5 & padj < 0.05 ~ "Down",
  log2FoldChange > 0.5 & padj < 0.05 ~ "Up",
  .default = "No Change"
)) %>% dplyr::mutate(labels = case_when(
  padj < 0.05 ~ gene,
  .default = ""
))

df_sifoxm1 <- df_sifoxm1 %>% na.omit()


summary <- df_sifoxm1 %>% group_by(direction) %>% summarise(n = n())

##Plot volcano, and change the number of up/down genes

ggplot(df_sifoxm1, aes(x = log2FoldChange, y = log10padj)) +
  ggrastr::rasterise(geom_point(aes(colour = direction), 
                                shape = 16, alpha = 0.5), dpi = 300) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c(
    "Down" = "blue",
    "No Change" = "gray",
    "Up" = "red"
  )) +
  ylim(0, 200) +
  xlim(-6.5, 6.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed") +
  labs(title = "Volcano Plot - 8y FOXM1 RNAi vs Mock",
       y = "-log10(padjust)", x = "Log2FC") +
  annotate("text", x = -6, y = 150, label = paste0("Down: ",summary[1,]$n), 
           hjust = "left") +
  annotate("text", x = -6, y = 140, label = paste0("Up: ", summary[3,]$n), 
           hjust = "left") +
  annotate("text", x = -6, y = 130, 
           label = paste0("No Change: ", summary[2,]$n), 
           hjust = "left") +
  annotate("rect", xmin=-6.2, xmax= -1, ymin=120 , 
           ymax=160, alpha=0.2, color="gray10", fill="white")



ggsave(
  "8y_sifoxm1_vs_mock_volcano_plot.pdf",
  path = "Plots",
  width = 8,
  height = 6,
  unit = "in"
)


## GSEA BPs

sifoxm1_ranked <- df_sifoxm1 %>% 
  dplyr::mutate(sign = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(log10pval = -log10(pvalue)) %>% 
  dplyr::mutate(metric = sign * log10pval) %>% 
  dplyr::arrange(desc(metric))

sifoxm1_ranked2 <- sifoxm1_ranked$metric
names(sifoxm1_ranked2) <- sifoxm1_ranked$gene
head(sifoxm1_ranked2, 50)


sifoxm1_hall <- GSEA(
  sifoxm1_ranked2,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-30,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  gson = NULL,
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea")


df_sifoxm1_hall <- as.data.frame(sifoxm1_hall)


df_sifoxm1_hall_plot <- df_sifoxm1_hall %>% 
  dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))

df_sifoxm1_hall_plot$Description <- gsub("HALLMARK_", "", 
                                         df_sifoxm1_hall_plot$Description)

ggplot(df_sifoxm1_hall_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-4, 4) + labs(y = "Hallmark Gene Sets",
                     title = "87y FOXM1 OE vs Vector",
                     subtitle = "Hallmark Collection") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")


ggsave(
  "8y_sifoxm1_vs_mock_hallmark.pdf",
  path = "Plots",
  width = 7,
  height = 8,
  unit = "in"
)



sifoxm1_gsea_bp <- gseGO(geneList = sifoxm1_ranked2,
                            OrgDb = org.Hs.eg.db,
                            ont  = "BP",
                            keyType = "SYMBOL",
                            pAdjustMethod = "BH",
                            minGSSize    = 30,
                            maxGSSize    = 700,
                            eps = 1e-300,
                            pvalueCutoff = 1,
                            verbose      = TRUE,
                            seed = TRUE,
                            by = "fgsea",
                            nPermSimple = 20000)


df_sifoxm1_gsea_bp <- as.data.frame(sifoxm1_gsea_bp)

nsig <- df_sifoxm1_gsea_bp %>% dplyr::filter(p.adjust < 0.25) %>% 
  dplyr::summarise(n = n()) %>% dplyr::pull(n)

sifoxm1_gsea_bp_005 <- sifoxm1_gsea_bp %>% 
  clusterProfiler::filter(p.adjust < 0.05)


df_sifoxm1_gsea_bp_plot <- df_sifoxm1_gsea_bp %>% 
  dplyr::mutate(log10padj = -log10(p.adjust)) %>% 
  dplyr::slice_max(log10padj, n = 30) %>% 
  dplyr::mutate(NES_sense = case_when(
    NES > 0 ~ "Activated",
    NES < 0 ~ "Supressed"
  ))


ggplot(df_sifoxm1_gsea_bp_plot, 
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(fill = NES_sense, size = p.adjust), shape = 21) + 
  theme_bw(base_size = 13) +
  xlim(-3, 3) + labs(y = "Gene Sets") +
  scale_size_continuous(range = c(7, 3), transform = "log10") +
  scale_fill_manual(values = c("Activated" = "red", "Supressed" = "blue"), 
                    name = "Direction")



network_data <- enrichmentNetwork(sifoxm1_gsea_bp_005@result, 
                                  verbose = T, 
                                  drawEllipses = TRUE, 
                                  fontSize = 2.5, 
                                  colorType = "nes",
                                  clustNameMethod = "pagerank",
                                  minClusterSize = 5,
                                  outerCutoff = 0.3) 

network_data + scale_color_gradient2(low = "blue", 
                                     mid = "white", 
                                     high = "red", midpoint = 0)

ggsave(
  "network_bps_8y_sifoxm1.pdf",
  path = "Plots",
  width = 7,
  height = 7,
  unit = "in"
)


## Merge siFOXM1 + All Other Data ----------------------------------------------





## Correlations FOXM1 OE -------------------------------------------------------

library(corrplot)
library(datasets)
library(heatmaply)
library(RColorBrewer)
library(gplots)
library(Hmisc)

oskm_genes <- partial_rep_vs_fibro$gene
aging_genes <- df_aging$gene_symbol
dndk_genes <- df_dndk$gene

common_genes <- intersect(oskm_genes, aging_genes)
common_genes <- intersect(common_genes, dndk_genes)

oskm_corr <- partial_rep_vs_fibro %>% 
  dplyr::select(gene, avg_log2FC) %>% 
  dplyr::rename(gene_symbol = gene, log2fc = avg_log2FC) %>% 
  dplyr::mutate(comparison = "oskm") %>% 
  dplyr::filter(gene_symbol %in% common_genes)

aging_corr <- df_aging %>% 
  dplyr::select(gene_symbol, log2FoldChange) %>% 
  dplyr::rename(gene_symbol = gene_symbol, log2fc = log2FoldChange) %>% 
  dplyr::mutate(comparison = "aging") %>% 
  dplyr::filter(gene_symbol %in% common_genes)


dndk_corr <- df_dndk %>% 
  dplyr::select(gene, log2FoldChange) %>% 
  dplyr::rename(gene_symbol = gene, log2fc = log2FoldChange) %>% 
  dplyr::mutate(comparison = "foxm1_oe") %>% 
  dplyr::filter(gene_symbol %in% common_genes)

corr_data <- rbind(oskm_corr, aging_corr, dndk_corr) %>% 
  pivot_wider(names_from = comparison, values_from = log2fc) %>% 
  as.data.frame()
row.names(corr_data) <- corr_data$gene_symbol
corr_data <- corr_data %>% dplyr::select(!gene_symbol)

res2 <- rcorr(as.matrix(corr_data))
res2
res2$r
res2$P


cor_matrix <- cor(corr_data, method = "spearman")
cor_matrix
corrplot(cor_matrix,
         method = "circle",
         type = "full",
         tl.cex = 0.8,
         is.corr = T,
         tl.col = "black",
         col = colorRampPalette(c("deepskyblue2", "white", "red"))(100),
         tl.srt = 45,
         hclust.method = c("complete"),
         addCoef.col = "black",
         addrect = 2,
         title = "Spearman Correlation Heatmap")



my_colors <- c("red", "blue")

col<- colorRampPalette(c("deepskyblue2", "white", "red"))(80)

heatmap.2(cor_matrix,  
          scale = "none",
          dendrogram = "column",
          distfun = dist,
          density.info = "none",
          trace = c("none"),
          col = col,
          colsep = c(1,2),
          rowsep = c(1,2),
          margins = c(8,8),
          cexRow = 1,
          cexCol = 1,
          revC = F,
          key = TRUE, 
          colRow = NULL,
          notecol = "black",
          hline =median(breaks),
          na.color=par("bg"),
          key.title = "Correlation",
          key.xlab = "Value")

## Correlation Log2FC


common_genes <- intersect(oskm_genes, dndk_genes)

partial_rep_vs_fibro_2 <- partial_rep_vs_fibro[common_genes,]

df_dndk_2 <- df_dndk[common_genes,]

row.names(df_dndk_2) == row.names(partial_rep_vs_fibro_2)

corr_rep_dndk <- data.frame(
  "genes" = df_dndk_2$gene,
  "FOXM1_OE" = df_dndk_2$log2FoldChange,
  "Partial_Reprogramming" = partial_rep_vs_fibro_2$avg_log2FC)

result_spearman = cor.test(corr_rep_dndk$FOXM1_OE, 
                           corr_rep_dndk$Partial_Reprogramming, 
                           method = "spearman")
result_spearman

library(ggpubr)

ggplot(corr_rep_dndk, aes(x = FOXM1_OE, y = Partial_Reprogramming)) +
  geom_point() +
  geom_smooth(method = "lm")




ggscatter(corr_rep_dndk, x = "FOXM1_OE", y = "Partial_Reprogramming", 
          add = "reg.line", conf.int = TRUE, shape = 16, alpha = 0.3,
          size = 1,
          cor.coef = TRUE,
          font.label = c(8, "bold", "black"),
          cor.coeff.args = list(method = "spearman", 
                                label.x = 1, label.sep = "\n"),
          xlab = "Log2FC FOXM1_OE", ylab = "Log2FC Partial_Reprogramming",
          add.params = list(color = "green", fill = "gray")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Log2FC Correlation (Spearman)")


ggsave(
  "dndk_reprog_single_gene_correlation.pdf",
  path = "Plots",
  width = 7,
  height = 7,
  unit = "in"
)

