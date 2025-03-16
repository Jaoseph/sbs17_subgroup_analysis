# ----- sysSVM2 output GSEA_enrichR

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
library(dplyr)
library(enrichR)

# Load driver predictions
predictions <- readRDS("test_sysSVM2/drivers_toppedUp.rds")
predictions <- subset(predictions, select = c(sample, canonical_driver, symbol))

# Load Sig17 annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample, Phenotype_Assigned))
names(occams_results_df) <- c("sample", "Phenotype_Assigned")

driver_df <- left_join(predictions, occams_results_df, by = "sample")
driver_df <- driver_df[, c(4,1,3,2)]

# Generate Matrix
matrix <- array(FALSE,c(length(unique(driver_df$symbol)),length(unique(driver_df$Phenotype_Assigned))))
colnames(matrix) <- unique(driver_df$Phenotype_Assigned)
rownames(matrix) <- unique(driver_df$symbol)

for(i in 1:nrow(driver_df)){
  matrix[driver_df[i,]$symbol,driver_df[i,]$Phenotype_Assigned] <- TRUE
}
matrix <- as.data.frame(matrix)

matrix$n_groups <- rowSums(matrix)
matrix$category <- apply(matrix[,c(1:4)], 1, function(x) paste(colnames(matrix[,c(1:4)])[x == TRUE], collapse = ":"))
matrix <- matrix %>%
  arrange(category)

BELike_drivers <- rownames(matrix)[matrix$category == "BarrettsLike.Sig17+"]
NaiveLike_drivers <- rownames(matrix)[matrix$category == "NaiveLike.Sig17+"]
TreatedLike_drivers <- rownames(matrix)[matrix$category == "TreatedLike.Sig17+"]
Negative_drivers <- rownames(matrix)[matrix$category == "Sig17-"]

# ------ Run EnrichR
dbs <- c("Reactome_2022")

# Run EnrichR --- BELike
BE_enriched <- enrichr(BELike_drivers, dbs)
BE_enriched <- BE_enriched$Reactome_2022
BE_enriched$Term <- gsub(" R-.*", "", BE_enriched$Term)
BE_enriched <- BE_enriched[BE_enriched$P.value < 0.1,]
BE_enriched <- BE_enriched %>%
  arrange(P.value) %>%
  mutate(Phenotype = "BarrettsLike.Sig17+")
#BE_enriched <- BE_enriched[c(1:10),]

pdf('Analysis/Group_Specific_Enrichment/BELike_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(BE_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Run EnrichR --- NaiveLike
Naive_enriched <- enrichr(NaiveLike_drivers, dbs)
Naive_enriched <- Naive_enriched$Reactome_2022
Naive_enriched$Term <- gsub(" R-.*", "", Naive_enriched$Term)
Naive_enriched <- Naive_enriched[Naive_enriched$P.value < 0.1,]
Naive_enriched <- Naive_enriched %>%
  arrange(P.value) %>%
  mutate(Phenotype = "NaiveLike.Sig17+")
#Naive_enriched <- Naive_enriched[c(1:10),]

pdf('Analysis/Group_Specific_Enrichment/NaiveLike_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(Naive_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Run EnrichR --- TreatedLike
Treated_enriched <- enrichr(TreatedLike_drivers, dbs)
Treated_enriched <- Treated_enriched$Reactome_2022
Treated_enriched$Term <- gsub(" R-.*", "", Treated_enriched$Term)
Treated_enriched <- Treated_enriched[Treated_enriched$P.value < 0.1,]
Treated_enriched <- Treated_enriched %>%
  arrange(P.value) %>%
  mutate(Phenotype = "TreatedLike.Sig17+")
#Treated_enriched <- Treated_enriched[c(1:10),]

pdf('Analysis/Group_Specific_Enrichment/TreatedLike_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(Treated_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
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

pdf('Analysis/Group_Specific_Enrichment/Negative_GSEAenrichR_Reactome.pdf', w = 8, h = 6)
plotEnrich(Negative_enriched, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()


# Generate Supplementary table
enrichR_output <- rbind(BE_enriched, Naive_enriched,Treated_enriched, Negative_enriched)
enrichR_output <- subset(enrichR_output, select = c("Phenotype", "Term", "Overlap", "P.value", "Genes"))

# Save
library(openxlsx)
write.xlsx(enrichR_output, "Analysis/Group_Specific_Enrichment/Reactome_Enrichment.xlsx", sheetName = "Table S")


# Non-overlapping MutSigDB

enriched_list <- list(NaiveLike = Naive_enriched, TreatedLike = Treated_enriched,
                      BarrettsLike = BE_enriched, Negative = Negative_enriched)


calculate_overlap <- function(group1, group2, term_column = "Term") {
  group1$Overlapping <- group1[[term_column]] %in% group2[[term_column]]
  group2$Overlapping <- group2[[term_column]] %in% group1[[term_column]]
  return(list(group1 = group1, group2 = group2))
}

pairwise_comparisons <- list(
  Naive_vs_Treated = calculate_overlap(enriched_list$NaiveLike, enriched_list$TreatedLike),
  Barretts_vs_Naive = calculate_overlap(enriched_list$BarrettsLike, enriched_list$NaiveLike),
  Naive_vs_Negative = calculate_overlap(enriched_list$NaiveLike, enriched_list$Negative),
  Barretts_vs_Negative = calculate_overlap(enriched_list$BarrettsLike, enriched_list$Negative))

nonoverlapping <- do.call(rbind, lapply(pairwise_comparisons, function(x) {
  rbind(transform(x$group1, Comparison = paste("Group1_vs_Group2")),
        transform(x$group2, Comparison = paste("Group2_vs_Group1")))
}))

nonoverlapping$Comparison <- rownames(nonoverlapping)
nonoverlapping$Comparison <- gsub("\\.\\d+$", "", nonoverlapping$Comparison)

# ----- Plot
library(tidyr)
nonoverlapping <- nonoverlapping %>%
  drop_na()
write.xlsx(nonoverlapping, "Analysis/Group_Specific_Enrichment/non_overlapping/Reactome/Reactome_Enrichment_nonoverlapping.xlsx", sheetName = "Table S")

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/Analysis/Group_Specific_Enrichment/")

# Only in Naive when comparing to Treated
NaiveNT <- subset(nonoverlapping, 
                  Phenotype == "NaiveLike.Sig17+" & Overlapping == FALSE & Comparison == "Naive_vs_Treated")
pdf('non_overlapping/Reactome/NaiveNT_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(NaiveNT, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Treated when comparing to Naive
TreatedNT <- subset(nonoverlapping, 
                    Phenotype == "TreatedLike.Sig17+" & Overlapping == FALSE & Comparison == "Naive_vs_Treated")
pdf('non_overlapping/Reactome/TreatedNT_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(TreatedNT, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Barretts when comparing to Naive
BarrettsBN <- subset(nonoverlapping, 
                     Phenotype == "BarrettsLike.Sig17+" & Overlapping == FALSE & Comparison == "Barretts_vs_Naive")
pdf('non_overlapping/Reactome/BarrettsBN_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(BarrettsBN, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Naive when comparing to Barretts
NaiveBN <- subset(nonoverlapping, 
                  Phenotype == "NaiveLike.Sig17+" & Overlapping == FALSE & Comparison == "Barretts_vs_Naive")
pdf('non_overlapping/Reactome/NaiveBN_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(NaiveBN, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Naive when comparing to Negative
NaiveNN <- subset(nonoverlapping, 
                  Phenotype == "NaiveLike.Sig17+" & Overlapping == FALSE & Comparison == "Naive_vs_Negative")
pdf('non_overlapping/Reactome/NaiveNN_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(NaiveNN, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Negative when comparing to Naive (optional)
NegativeNN <- subset(nonoverlapping, 
                     Phenotype == "Sig17-" & Overlapping == FALSE & Comparison == "Naive_vs_Negative")
pdf('non_overlapping/Reactome/NegativeNN_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(NegativeNN, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Barretts when comparing to Negative
BarrettsBNeg <- subset(nonoverlapping, 
                       Phenotype == "BarrettsLike.Sig17+" & Overlapping == FALSE & Comparison == "Barretts_vs_Negative")
pdf('non_overlapping/Reactome/BarrettsBNeg_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(BarrettsBNeg, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

# Only in Negative when comparing to Barretts (optional)
NegativeBNeg <- subset(nonoverlapping, 
                       Phenotype == "Sig17-" & Overlapping == FALSE & Comparison == "Barretts_vs_Negative")
pdf('non_overlapping/Reactome/NegativeBNeg_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
plotEnrich(NegativeBNeg, showTerms = 10, orderBy = "P.value", numChar = 70)
dev.off()

library(ggvenn)
library(tidyverse)
library(wesanderson)

venn_df <- split(nonoverlapping$Term, nonoverlapping$Phenotype)
venn_df <- venn_df[c("BarrettsLike.Sig17+", "NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-")]
#pdf('non_overlapping/Reactome/element_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
#ggvenn(venn_df, 
#       fill_color = c(wes_palette("AsteroidCity3")[3], wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[1], wes_palette("AsteroidCity3")[2]),
#       show_elements = TRUE, stroke_size = 0.5, set_name_size = 4, text_size = 2.3, label_sep = "\n") 
#dev.off()

#pdf('non_overlapping/Reactome/element_Reactome_Nonoverlapping.pdf', w = 8, h = 6)
#ggvenn(venn_df, 
#       fill_color = c(wes_palette("AsteroidCity3")[3], wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[1], wes_palette("AsteroidCity3")[2]),
#       show_elements = FALSE, stroke_size = 0.5, set_name_size = 4, text_size = 2.3) 


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
venn_df <- lapply(venn_df, wrap_terms, max_length = 20)

pdf('non_overlapping/Reactome/element_Reactome_Nonoverlapping_number.pdf', w = 10, h = 10)
ggvenn(venn_df, 
       fill_color = c(wes_palette("AsteroidCity3")[3], wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[1], wes_palette("AsteroidCity3")[2]),
       show_elements = TRUE, stroke_size = 0.5, set_name_size = 4, text_size = 4,label_sep = ",\n") 
dev.off()

