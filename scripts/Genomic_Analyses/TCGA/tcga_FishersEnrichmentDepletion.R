# ---- TCGA Fisher's Exact Test Depletion and Enrichment of Driver Genes

setwd('/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Genomic_Analysis')
library(rcompanion)
library(epitools)
library(dplyr)
library(ggpubr)

# Load Alteration Data for each classified group
load("tcga_mat_list.NaiveLike.Rdata") # List of matrices for NaiveLike samples
load("tcga_mat_list.Negative.Rdata")  # List of matrices for Sig17- samples

# Load EAC driver genes
drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

results <- data.frame()

alteration_type <- names(mat_list.n)

# Use alteration types from one group (assumes structure is identical across groups)
for (alteration in alteration_type) {
  for (driver in drivers) {
    if (driver %in% rownames(mat_list.p[[alteration]]) &&
        driver %in% rownames(mat_list.n[[alteration]])) {
      # Check if there are no mutations in both groups
      if (!any(mat_list.p[[alteration]][driver, ] == 1) &&
          !any(mat_list.n[[alteration]][driver, ] == 1)) {
        cat("Skipping driver", driver, "for alteration type", alteration, "- no mutations in both groups.\n")
      } else {
        # Generate contingency table
        groups <- list(
          NaiveLike = mat_list.p[[alteration]][driver, ],
          Negative = mat_list.n[[alteration]][driver, ]
        )
        contingency_table <- t(sapply(groups, function(x) c(n = sum(x == 0), m = sum(x == 1))))
        
        # Run Fisher's exact test
        fisher_test <- fisher.test(contingency_table)
        
        # Extract odds ratio and p-value
        odds <- as.data.frame(oddsratio.fisher(contingency_table)$measure)
        odds <- odds[odds$estimate != 1, ]
        odds$pval <- fisher_test$p.value
        
        # Add metadata
        odds$comparison <- "NaiveLike vs Negative"
        odds$reference <- "Negative"
        odds$driver <- driver
        odds$alteration <- alteration
        odds <- odds[, c("alteration", "driver", "comparison", "estimate", "lower", "upper", "pval")]
        
        results <- rbind(results, odds)
      }
    } else {
      cat("Driver", driver, "not found in alteration type", alteration, "\n")
    }
  }
}


# Adjust for multiple testing
results$padj <- p.adjust(results$pval, method = "BH")

# Save results
results <- results %>% arrange(padj, alteration)
write.csv(results, file = "EnrichmentDepletionComparison_TCGA.csv", row.names = FALSE)

# ----- Plot Significant Results ----- (No significant genes in TCGA)

# Filter results for significant findings
results <- results[results$padj < 0.1 & results$estimate != 0, ]
results$l2fc <- log2(results$estimate)
results$log2upper <- log2(results$upper)
results$log2lower <- log2(results$lower)
results$point_size <- abs(results$l2fc)
results$alteration <- factor(results$alteration, levels = c("SNV", "Indel", "AMP"))

pdf("EnrichmentDepletionComparison_TCGA.pdf")
ggplot(results, aes(x = l2fc, y = driver, xmin = log2lower, xmax = log2upper, color = alteration)) +
  geom_errorbarh(height = 0.2, color = "black", alpha = 0.7) +
  geom_point(size = results$point_size) +
  facet_wrap(~comparison, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "tomato") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "gray90", color = "black")
  ) +
  xlab("Log(Odds Ratio)") +
  scale_color_manual(name = "Alteration", values = c("#445E93", "#7EB2DD", "#F93943")) +
  scale_x_continuous(breaks = c(-7.5, -5, -2.5, 0, 2.5, 5))
dev.off()
