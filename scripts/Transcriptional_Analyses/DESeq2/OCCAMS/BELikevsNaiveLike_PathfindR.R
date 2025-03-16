# ---- Pathway Enrichment analysis BELike vs NaiveLike separately
library(pathfindR)
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/DESeq2/OCCAMS/")

# Load BE vs Negative significant DEGs
BENeg <- read.table(file ="BarrettsvsNegative_DESeq2_output.tsv")
BENeg <- BENeg[,c("Gene", "log2FoldChange", "padj")] ; rownames(BENeg) <- NULL

# Load Naive vs Negative significant DEGs
NaiveNeg <- read.table(file ="NaivevsNegative_DESeq2_output.tsv")
NaiveNeg <- NaiveNeg[,c("Gene", "log2FoldChange", "padj")] ; rownames(NaiveNeg) <- NULL

# Run PathfindR
BENeg.pathfindR <- run_pathfindR(BENeg, gene_sets = "Reactome",
                           min_gset_size = 5,
                           p_val_threshold = 0.1,
                           adj_method = "fdr",
                           enrichment_threshold = 0.05)
NaiveNeg.pathfindR <- run_pathfindR(NaiveNeg, gene_sets = "Reactome",
                                 min_gset_size = 5,
                                 p_val_threshold = 0.1,
                                 adj_method = "fdr",
                                 enrichment_threshold = 0.05)

# Clustering Enriched Terms
BENeg.pathfindR <- cluster_enriched_terms(BENeg.pathfindR, plot_dend = FALSE, plot_clusters_graph = FALSE)
BENeg.pathfindR <- BENeg.pathfindR[BENeg.pathfindR$Status == "Representative",]

NaiveNeg.pathfindR <- cluster_enriched_terms(NaiveNeg.pathfindR, plot_dend = FALSE, plot_clusters_graph = FALSE)
NaiveNeg.pathfindR <- NaiveNeg.pathfindR[NaiveNeg.pathfindR$Status == "Representative",]

# Save
library(openxlsx)
write.xlsx(BENeg.pathfindR, "PathfindR/BarrettsvsNegative_GSEApathfindR_Reactome.xlsx", rowNames = FALSE)
write.xlsx(NaiveNeg.pathfindR, "PathfindR/NaivesvsNegative_GSEApathfindR_Reactome.xlsx", rowNames = FALSE)

# Plot

BENeg.pathfindR <- BENeg.pathfindR[1:10,]
BENeg.pathfindR$Term_Description <- gsub("Regulation of Insulin-like Growth Factor \\(IGF\\) transport and uptake by Insulin-like Growth Factor Binding Proteins \\(IGFBPs\\)", 
                              "Regulation of IGF transport and uptake by IGF Binding Proteins", 
                              BENeg.pathfindR$Term_Description)


NaiveNeg.pathfindR <- NaiveNeg.pathfindR[1:10,]


pdf('PathfindR/BELikevsNegative_GSEApathfindR_Reactome.pdf')
enrichment_chart(
  result_df = BENeg.pathfindR,
  top_terms = 10)
dev.off()

pdf('PathfindR/NaiveLikevsNegative_GSEApathfindR_Reactome.pdf')
enrichment_chart(
  result_df = NaiveNeg.pathfindR,
  top_terms = 10)
dev.off()


# Test

load("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis/Processed_ExprMatrix.Rdata")
expr.matrix <- expr.matrix[expr.matrix$status %in% c("NaiveLike.Sig17+", "Sig17-"),]
Naive <-rownames(expr.matrix)[expr.matrix$status == "NaiveLike.Sig17+"]
Negative <- rownames(expr.matrix)[expr.matrix$status == "Sig17-"]
expr.matrix <- t(expr.matrix[,-1])

NaiveNeg.pathfindR <- NaiveNeg.pathfindR[1:10,]

pdf('PathfindR/NaiveLikevsNegative_score_matrix.pdf')
score_matrix <- score_terms(
  enrichment_table = NaiveNeg.pathfindR,
  exp_mat = expr.matrix,
  cases = Naive,
  use_description = TRUE, # default FALSE
  label_samples = FALSE, # default = TRUE
  case_title = "NaiveLike.Sig17+", # default = "Case"
  control_title = "Sig17-" # default = "Control"
  ,low = "#f7797d", # default = "green"
  ,mid = "#FFFFFF", # default = "black"
  high = "#1f4037" # default = "red"
)
dev.off()

# Test

load("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis/Processed_ExprMatrix.Rdata")
expr.matrix <- expr.matrix[expr.matrix$status %in% c("NaiveLike.Sig17+", "Sig17-"),]
Naive <-rownames(expr.matrix)[expr.matrix$status == "NaiveLike.Sig17+"]
Negative <- rownames(expr.matrix)[expr.matrix$status == "Sig17-"]
expr.matrix <- t(expr.matrix[,-1])

NaiveNeg.pathfindR <- NaiveNeg.pathfindR[1:10,]

pdf('PathfindR/NaiveLikevsNegative_score_matrix.pdf')
score_matrix <- score_terms(
  enrichment_table = NaiveNeg.pathfindR,
  exp_mat = expr.matrix,
  cases = Naive,
  use_description = TRUE, # default FALSE
  label_samples = FALSE, # default = TRUE
  case_title = "NaiveLike.Sig17+", # default = "Case"
  control_title = "Sig17-" # default = "Control"
  ,low = "#f7797d", # default = "green"
  ,mid = "#FFFFFF", # default = "black"
  high = "#1f4037" # default = "red"
)
dev.off()

# Test

load("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis/Processed_ExprMatrix.Rdata")
expr.matrix <- expr.matrix[expr.matrix$status %in% c("BarrettsLike.Sig17+"),]
BELike <-rownames(expr.matrix)[expr.matrix$status == "BarrettsLike.Sig17+"]
Negative <- rownames(expr.matrix)[expr.matrix$status == "Sig17-"]
expr.matrix <- t(expr.matrix[,-1])

NaiveNeg.pathfindR <- NaiveNeg.pathfindR[1:10,]

pdf('PathfindR/BELikevsNegativevsNegative_score_matrix.pdf')
score_matrix <- score_terms(
  enrichment_table = NaiveNeg.pathfindR,
  exp_mat = expr.matrix,
  cases = BELike,
  use_description = TRUE, # default FALSE
  label_samples = FALSE, # default = TRUE
  case_title = "BarrettsLike.Sig17+", # default = "Case"
  control_title = "Sig17-" # default = "Control"
  ,low = "#663333", # default = "green"
  ,mid = "#FFFFFF", # default = "black"
  high = "#222255" # default = "red"
)
dev.off()







