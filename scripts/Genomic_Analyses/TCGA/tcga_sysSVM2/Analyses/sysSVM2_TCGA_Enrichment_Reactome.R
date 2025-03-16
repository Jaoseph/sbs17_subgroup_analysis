# ----- sysSVM2 output GSEA_enrichR - tcga
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/Analysis")
library(dplyr)
library(devtools)
#install_github("wjawaid/enrichR")
library(enrichR)

# Load annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/sysSVM2_tcga_input.Rdata")

# A list of predicted drivers in each sample
driver_df <- readRDS("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/drivers_toppedUp.rds")
driver_df <- left_join(driver_df, molecular_data, by = c("sample", "entrez"))

# Load Sig17 annotations for tcga smaples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga$sample <- rownames(ann_tcga)
ann_tcga <- subset(ann_tcga, select = c(sample, Phenotype_Assigned))
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",] #remove BElike sample.

driver_df <- left_join(driver_df, ann_tcga, by = "sample")
driver_df <- driver_df[driver_df$sample != "TCGA-L5-A4OT",] #remove BElike sample.

# Generate Matrix
matrix <- array(FALSE,c(length(unique(driver_df$symbol)),length(unique(driver_df$Phenotype_Assigned))))
colnames(matrix) <- unique(driver_df$Phenotype_Assigned)
rownames(matrix) <- unique(driver_df$symbol)

for(i in 1:nrow(driver_df)){
  matrix[driver_df[i,]$symbol,driver_df[i,]$Phenotype_Assigned] <- TRUE
}
matrix <- as.data.frame(matrix)

matrix$n_groups <- rowSums(matrix)
matrix$category <- apply(matrix[,c(1:2)], 1, function(x) paste(colnames(matrix[,c(1:2)])[x == TRUE], collapse = ":"))
matrix <- matrix %>%
  arrange(category)

NaiveLike_drivers <- rownames(matrix)[matrix$category == "NaiveLike.Sig17+"]
Negative_drivers <- rownames(matrix)[matrix$category == "Sig17-"]


# ------ Run EnrichR
dbs <- c("Reactome_2022")

# Run EnrichR --- NaiveLike
Naive_enriched <- enrichr(NaiveLike_drivers, dbs)
Naive_enriched <- Naive_enriched$Reactome_2022
Naive_enriched$Term <- gsub(" R-.*", "", Naive_enriched$Term)
Naive_enriched <- Naive_enriched[Naive_enriched$P.value < 0.1,]
Naive_enriched <- Naive_enriched %>%
  arrange(P.value) %>%
  mutate(Phenotype = "NaiveLike.Sig17+")
#Naive_enriched <- Naive_enriched[c(1:10),]

pdf('Group_Specific_Enrichment/NaiveLike_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(Naive_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Run EnrichR --- Negative
Negative_enriched <- enrichr(Negative_drivers, dbs)
Negative_enriched <- Negative_enriched$Reactome_2022
Negative_enriched$Term <- gsub(" R-.*", "", Negative_enriched$Term)
Negative_enriched <- Negative_enriched[Negative_enriched$P.value < 0.1,]
Negative_enriched <- Negative_enriched %>%
  arrange(P.value) %>%
  mutate(Phenotype = "Sig17-")
#Negative_enriched <- Negative_enriched[c(1:10),]

pdf('Group_Specific_Enrichment/Negative_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(Negative_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()


# Generate Supplementary table
Naive_enriched$Overlapping <- ifelse(Naive_enriched$Term %in% Negative_enriched$Term, TRUE, FALSE)
Negative_enriched$Overlapping <- ifelse(Negative_enriched$Term %in% Naive_enriched$Term, TRUE, FALSE)

enrichR_output <- rbind(Naive_enriched, Negative_enriched)
enrichR_output <- subset(enrichR_output, select = c("Phenotype", "Term", "Overlap", "P.value", "Genes", "Overlapping"))

# Save
library(openxlsx)
write.xlsx(enrichR_output, "Group_Specific_Enrichment/Reactome_Enrichment.xlsx", sheetName = "STab 16. TCGA Drivers EnrichR")

# Plot Naive Non-overlapping Reactome Pathways
Naive_nonoverlapping <- enrichR_output[enrichR_output$Phenotype == "NaiveLike.Sig17+" & enrichR_output$Overlapping == FALSE,]

pdf('Group_Specific_Enrichment/NaiveLike_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(Naive_nonoverlapping, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()



# Test

library(ggvenn)
library(tidyverse)
library(wesanderson)

enriched_list <- list(NaiveLike = Naive_enriched, Negative = Negative_enriched)


calculate_overlap <- function(group1, group2, term_column = "Term") {
  group1$Overlapping <- group1[[term_column]] %in% group2[[term_column]]
  group2$Overlapping <- group2[[term_column]] %in% group1[[term_column]]
  return(list(group1 = group1, group2 = group2))
}

pairwise_comparisons <- list(Naive_vs_Negative = calculate_overlap(enriched_list$NaiveLike, enriched_list$Negative))

nonoverlapping <- do.call(rbind, lapply(pairwise_comparisons, function(x) {
  rbind(transform(x$group1, Comparison = paste("Group1_vs_Group2")),
        transform(x$group2, Comparison = paste("Group2_vs_Group1")))
}))

nonoverlapping$Comparison <- rownames(nonoverlapping)
nonoverlapping$Comparison <- gsub("\\.\\d+$", "", nonoverlapping$Comparison)

venn_df <- split(nonoverlapping$Term, nonoverlapping$Phenotype)
venn_df <- venn_df[c("NaiveLike.Sig17+","Sig17-")]

wrap_terms <- function(terms, max_length = 30) {
  sapply(terms, function(term) {
    if (nchar(term) > max_length) {
      # Split term by words or characters and add newlines
      paste(strwrap(term, width = max_length), collapse = "\n")
    } else {
      term
    }
  })
}
venn_df <- lapply(venn_df, wrap_terms, max_length = 30)

pdf('Group_Specific_Enrichment/element_Reactome_Nonoverlapping.pdf', w = 10, h = 10)
ggvenn(venn_df,
       fill_color = c(wes_palette("AsteroidCity3")[3], wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[1], wes_palette("AsteroidCity3")[2]),
       show_elements = TRUE, stroke_size = 0.5, set_name_size = 5, text_size = 2.3, label_sep = ",\n") 
dev.off()

