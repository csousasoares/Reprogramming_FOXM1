## Load Packages and Set Seed --------------------------------------------------

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(msigdbr)
library(enrichplot)


## Load 2022 In Vivo Reprog Belmonte Data --------------------------------------

in_vivo_7m_belmonte_counts <- read.csv("Input Files\\2022_in_vivo_reprog_belmonte\\in_vivo_reprog_belmonte_seven_month_counts.csv",
                                       header = T, sep = ";", row.names = 1)

in_vivo_7m_belmonte_metadata <- read.csv("Input Files\\2022_in_vivo_reprog_belmonte\\in_vivo_reprog_belmonte_seven_month_metadata.csv",
                                         sep = ";", header = T) ## Curated 

in_vivo_7m_belmonte_metadata_2 <- read.csv("Input Files\\2022_in_vivo_reprog_belmonte\\in_vivo_reprog_belmonte_seven_month_metadata_2.csv",
                                           sep = ";", header = T) ## From GEO

## Sanity check if everything matches:

in_vivo_7m_belmonte_metadata %>% 
  dplyr::arrange(Sample) %>% pull(Sample) == in_vivo_7m_belmonte_metadata_2 %>% 
  dplyr::arrange(Sample) %>% pull(Sample) ## Should be TRUE


## Get the full metadata:

in_vivo_7m_belmonte_full_meta <- in_vivo_7m_belmonte_metadata %>% 
  dplyr::left_join(in_vivo_7m_belmonte_metadata_2, by = "Sample") ## Join meta

row.names(in_vivo_7m_belmonte_full_meta) <- in_vivo_7m_belmonte_full_meta$Sample

in_vivo_7m_belmonte_full_meta <- in_vivo_7m_belmonte_full_meta %>% 
  dplyr::select(!Sample) ## Put in correct format, samples as rownames


liver_samples <- row.names(in_vivo_7m_belmonte_full_meta) %>% 
  str_detect("Liver")

liver_samples <- row.names(in_vivo_7m_belmonte_full_meta[liver_samples,])

liver_2022_counts <- in_vivo_7m_belmonte_counts[,liver_samples]
liver_2022_meta <- in_vivo_7m_belmonte_full_meta[liver_samples,]

colnames(liver_2022_counts) == row.names(liver_2022_meta) ## Should be True

dds_sex <- DESeqDataSetFromMatrix(countData = liver_2022_counts,
                                  colData = liver_2022_meta,
                                  design = ~ Sex + Condition_2)
dds_sex


smallestGroupSize <- 10
keep <- rowSums(counts(dds_sex) >=10) >= smallestGroupSize ##Recomended in vignette
dds_sex <- dds_sex[keep,]

dds_sex <- DESeq(dds_sex)

result_4f_liver <- results(dds_sex, contrast = c("Condition_2", "4F", "Control"), 
                        cooksCutoff=FALSE, independentFiltering=FALSE) 

df_4f_liver_sex <- as.data.frame(result_4f_liver)
df_4f_liver_sex$gene <- row.names(df_4f_liver_sex)

annots <- AnnotationDbi::select(org.Mm.eg.db, keys = df_4f_liver_sex$gene,
                                columns="SYMBOL", keytype="ENSEMBL")

annots_summary <- annots %>% 
  na.omit() %>% 
  dplyr::group_by(ENSEMBL) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n == 1)

df_4f_liver_sex_filtered <- df_4f_liver_sex %>% 
  dplyr::filter(gene %in% annots_summary$ENSEMBL)

annots <- AnnotationDbi::select(org.Mm.eg.db, keys = df_4f_liver_sex_filtered$gene,
                                columns="SYMBOL", keytype="ENSEMBL")


df_4f_liver_sex_filtered$gene_symbol <- annots$SYMBOL

## Most significant is Serpina3c aka ENSMUSG00000066361, with log2fc ~ 7

liver_2022_counts["ENSMUSG00000066361",] ## It looks indeed highly up


## GSEA Hallmarks 4F Liver 7m --------------------------------------------------


## Create ranked list:

df_4f_liver_sex_filtered <- df_4f_liver_sex_filtered %>% 
  dplyr::mutate(direction = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(log10pval = -log10(pvalue)) %>% 
  dplyr::mutate(metric = direction * log10pval) %>% 
  dplyr::arrange(desc(metric))

df_4f_liver_sex_ranked <- df_4f_liver_sex_filtered$metric
names(df_4f_liver_sex_ranked) <- df_4f_liver_sex_filtered$gene_symbol
head(df_4f_liver_sex_ranked, 10)
tail(df_4f_liver_sex_ranked, 10)

## Create Hallmarks Term to Gene:

Hs <- msigdbr(
  species = "Mus musculus",
  collection = "H")
Hs

Hs_t2g = Hs %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
Hs_t2g


GSEA_4f_liver_H <- GSEA(df_4f_liver_sex_ranked,
                  exponent = 1,
                  minGSSize = 10,
                  maxGSSize = 1000,
                  eps = 1e-300,
                  pvalueCutoff = 1,
                  pAdjustMethod = "BH",
                  gson = NULL,
                  TERM2GENE = Hs_t2g,
                  verbose = T,
                  seed = T,
                  by = "fgsea")

results_GSEA_4f_liver_H <- as.data.frame(GSEA_4f_liver_H)

gseaplot2(GSEA_4f_liver_H, geneSetID = c(1,2,3),
          base_size = 12,
          color = "green",
          rel_heights = c(2,0.5,0.5),
          subplots = 1:3,
          pvalue_table = T,
          ES_geom = "line")


## Write session_info:

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

