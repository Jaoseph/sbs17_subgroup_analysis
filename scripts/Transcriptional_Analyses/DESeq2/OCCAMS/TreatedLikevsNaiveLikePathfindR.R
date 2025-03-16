# ---- Pathway Enrichment analysis (Treated vs Naive) for OCCAMS only
library(pathfindR)

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/DESeq2/OCCAMS/")

# Load Treated vs Naive significant DEGs
gene_list <- read.table(file ="TreatedvsNaive_DESeq2_output.tsv")
gene_list <- gene_list[,c("Gene", "log2FoldChange", "padj")] ; rownames(gene_list) <- NULL

# Run PathfindR
pathfindR <- run_pathfindR(gene_list, gene_sets = "Reactome",
                                 min_gset_size = 5,
                                p_val_threshold = 0.1,
                                 adj_method = "fdr",
                                 enrichment_threshold = 0.05)

# Clustering Enriched Terms
clustered_output_df <- cluster_enriched_terms(pathfindR, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_output_df <- clustered_output_df[clustered_output_df$Status == "Representative",]

# Save
pdf('PathfindR/TreatedvsNaive_GSEApathfindR_Reactome.pdf')
enrichment_chart(
  result_df = clustered_output_df,
  top_terms = 10)
dev.off()

library(openxlsx)
write.xlsx(clustered_output_df, "PathfindR/TreatedvsNaive_GSEApathfindR_Reactome.xlsx", rowNames = FALSE)

# Test

load("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis/Processed_ExprMatrix.Rdata")
expr.matrix <- expr.matrix[expr.matrix$status %in% c("NaiveLike.Sig17+", "TreatedLike.Sig17+"),]
Treated <- rownames(expr.matrix)[expr.matrix$status == "TreatedLike.Sig17+"]
Naive <-rownames(expr.matrix)[expr.matrix$status == "NaiveLike.Sig17+"]
expr.matrix <- t(expr.matrix[,-1])

clustered_output_df <- clustered_output_df[1:10,]

pdf('PathfindR/TreatedvsNaive_score_matrix.pdf')
score_matrix <- score_terms(
  enrichment_table = clustered_output_df,
  exp_mat = expr.matrix,
  cases = Treated,
  use_description = TRUE, # default FALSE
  label_samples = FALSE, # default = TRUE
  case_title = "TreatedLike.Sig17+", # default = "Case"
  control_title = "NaiveLike.Sig17+", # default = "Control"
  low = "#f7797d", # default = "green"
  mid = "#fffde4", # default = "black"
  high = "#1f4037" # default = "red"
)
dev.off()
