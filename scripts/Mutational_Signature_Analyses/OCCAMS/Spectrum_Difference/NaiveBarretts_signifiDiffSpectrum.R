##### Comparing mutational signature contribution patterns

setwd("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/Statistical_testing")
library(dplyr)
library(ggpubr)

# ----- Naive vs Barrett

## Load Mean Representative Mutational Spectra

#Naive
load("naive_mut_dists.Rdata")
load('/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Naive_IDnormalised_PhenotypeAnnotation_clust.Rdata')
naive.ann <- ann[ann$project_code == "ESAD" & ann$Phenotype == "Sig17+",]
naive.mut.dists <-  naive.mut.dists[rownames(naive.mut.dists) %in% rownames(naive.ann),]

#Barrett
load("barretts_mut_dists.Rdata")
load('/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata')
ann$Phenotype <- as.character(ann$Phenotype)
barretts.ann <- ann[ann$Phenotype == "Sig17+", ,drop = FALSE]
barretts.mut.dists <- barretts.mut.dists[rownames(barretts.mut.dists) %in% rownames(barretts.ann),]

rm(ann)

# For each trinucleotide and indel context, I want to run wilcoxon test, comparing the Naive and Barrett samples

wilcoxon_output <- data.frame(matrix(nrow = ncol(naive.mut.dists), ncol = 3))
colnames(wilcoxon_output) <- c("Context", "Wstat", "Pval")

for (i in seq_len(ncol(naive.mut.dists))) {
  
  current_context <- colnames(naive.mut.dists)[i]
  
  naive_contributions <- naive.mut.dists[,current_context]
  barretts_contributions <- barretts.mut.dists[,current_context]
  
  # Run statistical test
  result <- wilcox.test(naive_contributions, barretts_contributions)
  
  # Store test results
  wilcoxon_output$Context[i] <- current_context
  wilcoxon_output$Wstat[i] <- result$statistic
  wilcoxon_output$Pval[i] <- result$p.value
  
}

wilcoxon_output$padj <- p.adjust(wilcoxon_output$Pval, method = "BH")
wilcoxon_output <- wilcoxon_output %>%
  arrange(padj)

wilcoxon_output <- wilcoxon_output[wilcoxon_output$padj < 0.05,]

# Save
library(openxlsx)
write.xlsx(wilcoxon_output, "NaiveBarrett_signifiDiffSpectrum.xlsx", sheetName = "Table S", rowNames=TRUE)