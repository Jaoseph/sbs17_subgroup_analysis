# ----- Compare DEGs from BarrettsLike.Sig17+ vs Negative AND Naivelike.Sig17+ vs Negative 
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
BN <- read.table("DESeq2/OCCAMS/BarrettsvsNegative_DESeq2_output.tsv") 
NN <- read.table("DESeq2/OCCAMS/NaivevsNegative_DESeq2_output.tsv")
common_degs <- intersect(BN$Gene, NN$Gene)
write.table(common_degs, file = "NaiveLikevsBarrettsLike_commonDEGS.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


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
write.xlsx(collated, "NaiveLikevsBarrettsLike_commonDEGS.xlsx", rowNames = FALSE)

# ----- Plot

collated <- collated %>%
  mutate(signif = ifelse(BN_log2FoldChange > 0 & NN_log2FoldChange > 0, "Up-regulated", 
                         ifelse(BN_log2FoldChange < 0 & NN_log2FoldChange < 0, "Down-regulated", "Mixed")))

pdf("NaiveLikevsBarrettsLike_commonDEGS.pdf",w=6,h=6)
ggscatter(collated, y = "NN_log2FoldChange", x = "BN_log2FoldChange",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"), color = "signif", 
          label = "Gene", repel = TRUE, ylab = 'Log2 fold change NaiveLikevsNegative', xlab = "Log2 fold change BarrettsLikevsNegative") +  
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth=0.3)+
  geom_vline(xintercept=0, linetype="dashed", color = "black", linewidth=0.3)
dev.off()


