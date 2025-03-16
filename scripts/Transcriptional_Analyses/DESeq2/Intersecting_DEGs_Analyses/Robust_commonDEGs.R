# ----- consistent across subtypes(BELike and NaiveLike) AND across cohorts (TCGA and OCCAMS)
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis")

# Load Libraries
library(biomaRt)
library(data.table)
library(dplyr)
library(ggpubr)
library(tidyverse)


# Connect to BioMart database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Load Significant DEGS
NB <- read.delim("NaiveLikevsBarrettsLike_commonDEGS.txt", header=FALSE)$V1
OT <- read.delim("occams_tcga_commonDEGS.txt", header=FALSE)$V1
common_degs <- intersect(NB, OT)
write.table(common_degs, file = "Robust_commondegs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

BN <- read.table("DESeq2/OCCAMS/BarrettsvsNegative_DESeq2_output.tsv") 
NN <- read.table("DESeq2/OCCAMS/NaivevsNegative_DESeq2_output.tsv")


# Re-format BN DESeq2 ouput
BN <- BN %>%
  mutate(BN_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% common_degs) %>%
  dplyr::select(Gene, BN_log2FoldChange) %>%
  arrange(desc(BN_log2FoldChange)) ; rownames(BN) <- NULL

# Re-format NN DESeq2 ouput
NN <- NN %>%
  mutate(NN_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% common_degs) %>%
  dplyr::select(Gene, NN_log2FoldChange) %>%
  arrange(desc(NN_log2FoldChange)) ; rownames(NN) <- NULL

# Combine DESeq2 outputs

gene_info <- getBM(attributes = c("hgnc_symbol", "description"),
                   filters = "hgnc_symbol",
                   values = common_degs,
                   mart = ensembl)
gene_info$description <- gsub("\\[.*?\\]", "", gene_info$description)
colnames(gene_info) <- c("Gene", "Description")

collated <- merge(BN, NN, by = "Gene")
collated <- merge(collated, gene_info, by = "Gene")
collated <- collated %>%
  arrange(Gene)

library(openxlsx)
write.xlsx(collated, "Robust_commondegs.xlsx", rowNames = FALSE)

# ----- Plot

collated <- collated %>%
  mutate(signif = ifelse(BN_log2FoldChange > 0 & NN_log2FoldChange > 0, "Up-regulated", 
                         ifelse(BN_log2FoldChange < 0 & NN_log2FoldChange < 0, "Down-regulated", "Mixed")))



#custom label
custom_labels <- c(
  "Log2 fold change NaiveLikevsNegative" = "l2fc NaiveLike vs **Negative**",
  "Log2 fold change BarrettsLikevsNegative" = "Barrettslike vs **Negative**")

labeller = as_labeller(custom_labels)

pdf("Robust_commondegs.pdf",w=5,h=5)
ggscatter(collated, y = "NN_log2FoldChange", x = "BN_log2FoldChange",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"), color = "signif", 
          label = "Gene", repel = TRUE, ylab = 'l2fc NaiveLikevsNegative', xlab = "l2fc BarrettsLikevsNegative") +  
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", linewidth=0.3) +
  labs(color = '')
dev.off()

# --- TCGA vs OCCAMS

# Load Significant DEGS
NB <- read.delim("NaiveLikevsBarrettsLike_commonDEGS.txt", header=FALSE)$V1
OT <- read.delim("occams_tcga_commonDEGS.txt", header=FALSE)$V1
common_degs <- intersect(NB, OT)

occams <- read.table("DESeq2/OCCAMS/NaivevsNegative_DESeq2_output.tsv")
tcga <- read.table("DESeq2/TCGA/NaivevsNegative_DESeq2_output.tsv")


# Re-format OCCAMS DESeq2 ouput
occams <- occams %>%
  mutate(occams_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% common_degs) %>%
  dplyr::select(Gene, occams_log2FoldChange) %>%
  arrange(desc(occams_log2FoldChange)) ; rownames(occams) <- NULL

# Re-format TCGA DESeq2 ouput
tcga <- tcga %>%
  mutate(tcga_log2FoldChange = log2FoldChange) %>%
  filter(Gene %in% common_degs) %>%
  dplyr::select(Gene, tcga_log2FoldChange) %>%
  arrange(desc(tcga_log2FoldChange)) ; rownames(tcga) <- NULL

# Combine DESeq2 outputs

gene_info <- getBM(attributes = c("hgnc_symbol", "description"),
                   filters = "hgnc_symbol",
                   values = common_degs,
                   mart = ensembl)
gene_info$description <- gsub("\\[.*?\\]", "", gene_info$description)
colnames(gene_info) <- c("Gene", "Description")

collated <- merge(tcga, occams, by = "Gene")
collated <- merge(collated, gene_info, by = "Gene")
collated <- collated %>%
  arrange(Gene)

library(openxlsx)
write.xlsx(collated, "Robust_commondegs_TCGAOCCAMS.xlsx", rowNames = FALSE)

# ----- Plot

collated <- collated %>%
  mutate(signif = ifelse(occams_log2FoldChange > 0 & tcga_log2FoldChange > 0, "Up-regulated", 
                         ifelse(occams_log2FoldChange < 0 & tcga_log2FoldChange < 0, "Down-regulated", "Mixed")))

pdf("Robust_commondegs_TCGAOCCAMS.pdf",w=5,h=5)
ggscatter(collated, x = "tcga_log2FoldChange", y = "occams_log2FoldChange",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"), color = "signif", 
          label = "Gene", repel = TRUE, xlab = 'l2fc NaiveLike TCGA', ylab = "l2fc NaiveLike OCCAMS") +  
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", linewidth=0.3)+
  labs(color = '')
dev.off()
