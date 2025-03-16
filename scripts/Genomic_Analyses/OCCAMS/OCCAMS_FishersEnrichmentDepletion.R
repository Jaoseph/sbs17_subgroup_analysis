# ---- OCCAMS Fisher's Exact Test Depletion and Enrichment of Driver Genes

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/genomic_changes/Oncoprint")
library(rcompanion)
library(epitools)
library(dplyr)
library(ggpubr)

# Load Alteration Data for each classified group

load("mat_list.BELike.Rdata")
load("mat_list.NaiveLike.Rdata")
load("mat_list.TreatedLike.Rdata")
load("mat_list.Negative.Rdata")

# Load EAC driver genes
drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

results <- data.frame()

alteration_type <- names(mat_list.BELike)

for (alteration in alteration_type){
  for (driver in drivers) {
    #generate contingency table
    groups <- list(
      BElike = mat_list.BELike[[alteration]][driver,],
      Naivelike = mat_list.NaiveLike[[alteration]][driver,],
      Treatedlike = mat_list.TreatedLike[[alteration]][driver,],
      Negative = mat_list.Negative[[alteration]][driver,]
    )
    contingency_table <- t(sapply(groups, function(x) c(n = sum(x == 0),m = sum(x == 1))))
    
    #run pairwise-fishers test
    PT <- pairwiseNominalIndependence(contingency_table,
                                      fisher = TRUE,
                                      gtest = FALSE,
                                      chisq = FALSE,
                                      method = "none",
                                      simulate.p.value = TRUE)
    
    #remove BElike vs Treated and TreatedLike vs Negative (quite uninformative, only increases multiple testing burden)
    PT <- PT[!(PT$Comparison %in% c("BElike : Treatedlike", "Treatedlike : Negative")),] 
    
    for (comparison in PT$Comparison) {
      #extract the two groups with significance
      groups <- strsplit(comparison, " : ")[[1]]
      groupA <- groups[1]
      groupB <- groups[2]
      
      if (groupA == "BElike" & groupB == "Naivelike"){
        ref <- groupB
        non_ref <- groupA
      } else if (groupA == "BElike" & groupB == "Negative"){
        ref <- groupB 
        non_ref <- groupA
      } else if (groupA == "Naivelike" & groupB == "Treatedlike"){
        ref <- groupA
        non_ref <- groupB
      } else if (groupA == "Naivelike" & groupB == "Negative"){
        ref <- groupB
        non_ref <- groupA
      }
      
      #sub-set contigency table from original contingency table matrix
      submatrix <- contingency_table[c(ref, non_ref),]
      
      #calculate odds ratio
      odds <- as.data.frame(oddsratio.fisher(submatrix)$measure)
      odds <- odds[odds$estimate != 1, ] 
      odds$pval <- PT$p.Fisher[PT$Comparison == comparison]
      
      odds$comparison <- comparison
      odds$reference <- ref
      
      odds$driver <- driver
      odds <- odds[,c("driver", setdiff(names(odds), "driver"))]
      
      odds$alteration <- alteration
      odds <- odds[,c("alteration", setdiff(names(odds), "alteration"))]
      
      results <- rbind(results, odds)
      rownames(results) <- NULL
    }
  }
}

# Adjust for P-value
results$padj <- p.adjust(results$pval, method = "BH")

# Save
results <- results %>%
  arrange(padj, alteration)

results$comparison <- gsub(" : ", " vs ", results$comparison)

#write.csv(results, file = "EnrichmentDepletionComparison.csv", row.names = F)


# ----- Plot

results <- results[results$padj < 0.1,]
results <- results[results$estimate != 0,]
results$l2fc <- log2(results$estimate)
results$log2upper <- log2(results$upper)
results$log2lower <- log2(results$lower)
results$point_size <- abs(results$l2fc) 
results$alteration <- factor(results$alteration, levels = c("SNV", "Indel", "AMP"))

#custom label
custom_labels <- c(
  "BElike vs Naivelike" = "BElike vs **Naivelike**",
  "BElike vs Negative" = "BElike vs **Negative**",
  "Naivelike vs Negative" = "Naivelike vs **Negative**")

pdf("EnrichmentDepletionComparison.pdf", w = 6, h = 6)
ggplot(results, aes(x = l2fc, y = driver, xmin = log2lower, xmax = log2upper, color = alteration)) +
  geom_errorbarh(height = 0.2, color = "black", alpha = 0.7) +
  geom_point(size = results$point_size) +
  facet_wrap(~comparison, nrow = 1, labeller = as_labeller(custom_labels)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#CC3311") +
  theme_minimal() +
  theme(axis.title.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "gray90", color = "black"), strip.text = ggtext::element_markdown()) + xlab("Log(Odds Ratio)") +
  scale_color_manual(name = "Alteration", values = c("#FFAABB","#332288","#999933")) +
  scale_x_continuous(breaks = c(-7.5,-5,-2.5,0,2.5,5))
dev.off()

