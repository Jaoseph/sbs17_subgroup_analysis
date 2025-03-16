# ----- Find common genes differentially expressed in both OCCAMs and TCGA 
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis")

# Load Libraries
library(biomaRt)
library(data.table)
library(dplyr)
library(ggpubr)
library(tidyverse)


# Connect to BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Identify Common DEGS across cohorts
occams <- read.table("DESeq2/OCCAMS/NaivevsNegative_DESeq2_output.tsv")
tcga <- read.table("DESeq2/TCGA/NaivevsNegative_DESeq2_output.tsv")
common_cohort <- intersect(occams$Gene, tcga$Gene)

# Identify Common DEGS across BELike and NaiveLike in OCCAMS
NB <- read.delim("NaiveLikevsBarrettsLike_commonDEGS.txt", header=FALSE)$V1
OT <- read.delim("occams_tcga_commonDEGS.txt", header=FALSE)$V1
common_subgroup <- intersect(NB, OT)

# Find difference (i.e. genes specific to Naivelike, common across cohorts)
diff_degs <- setdiff(common_cohort, common_subgroup)


# Re-format OCCAMS DESeq2 ouput
occams <- occams %>%
  mutate(occams_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% diff_degs) %>%
  dplyr::select(Gene, occams_log2FoldChange) %>%
  arrange(desc(occams_log2FoldChange)) ; rownames(occams) <- NULL

# Re-format TCGA DESeq2 ouput
tcga <- tcga %>%
  mutate(tcga_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% diff_degs) %>%
  dplyr::select(Gene, tcga_log2FoldChange) %>%
  arrange(desc(tcga_log2FoldChange)) ; rownames(tcga) <- NULL

# Combine DESeq2 outputs

gene_info <- getBM(attributes = c("hgnc_symbol", "description"),
                   filters = "hgnc_symbol",
                   values = diff_degs,
                   mart = ensembl)
gene_info$description <- gsub("\\[.*?\\]", "", gene_info$description)
colnames(gene_info) <- c("Gene", "Description")

collated <- merge(tcga, occams, by = "Gene")
collated <- merge(collated, gene_info, by = "Gene")
collated <- collated %>%
  arrange(Gene)

library(openxlsx)
write.xlsx(collated, "occams_tcga_NaiveLikespecific_commonDEGS.xlsx", rowNames = FALSE)

# ----- Plot

collated <- collated %>%
  mutate(signif = ifelse(occams_log2FoldChange > 0 & tcga_log2FoldChange > 0, "Up-regulated", 
                         ifelse(occams_log2FoldChange < 0 & tcga_log2FoldChange < 0, "Down-regulated", "Mixed")))

pdf("occams_tcga_NaiveLikespecific_commonDEGS.pdf",w=6,h=6)
ggscatter(collated, x = "tcga_log2FoldChange", y = "occams_log2FoldChange",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"), color = "signif", 
          label = "Gene", repel = TRUE, xlab = 'Log2 fold change TCGA', ylab = "Log2 fold change OCCAMS") +  
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", linewidth=0.3)
dev.off()

