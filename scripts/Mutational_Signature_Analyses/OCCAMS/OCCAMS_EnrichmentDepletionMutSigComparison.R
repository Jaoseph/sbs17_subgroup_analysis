# ---- OCCAMS Fisher's Exact Test Depletion and Enrichment of Mutational Signatures

setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs")

# Load DeconstructSigs contribution and change into binary form
load("Primaries_deconstructSigs.Rdata")
primaries_sigs_complete <- ifelse(primaries_sigs_complete > 0.05, 1, 0)

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample ; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)
primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names ; primaries_sigs_complete <- primaries_sigs_complete[,-1]

primaries_sigs_complete <- subset(primaries_sigs_complete, select = -c(SBS17a, SBS17b))

# Loop to run pairwise fishers

results <- data.frame()

MutSigs <- colnames(primaries_sigs_complete[,-1])

for (sig in MutSigs) {
  groups <- list(
    BElike = primaries_sigs_complete %>%
      filter(Phenotype_Assigned == "BarrettsLike.Sig17+") %>%
      select(all_of(sig)),
    Naivelike = primaries_sigs_complete %>%
      filter(Phenotype_Assigned == "NaiveLike.Sig17+") %>%
      select(all_of(sig)),
    Treatedlike = primaries_sigs_complete %>%
      filter(Phenotype_Assigned == "TreatedLike.Sig17+") %>%
      select(all_of(sig)),
    Negative = primaries_sigs_complete %>%
      filter(Phenotype_Assigned == "Sig17-") %>%
      select(all_of(sig))
  )
  contingency_table <- t(sapply(groups, function(x) c(absent = sum(x == 0), present = sum(x == 1))))
  
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
    
    odds$Sig <- sig
    odds <- odds[,c("Sig", setdiff(names(odds), "Sig"))]

    results <- rbind(results, odds)
    rownames(results) <- NULL
  }
}

# Adjust for P-value
results$padj <- p.adjust(results$pval, method = "BH")

results <- results %>%
  arrange(padj)
results$comparison <- gsub(" : ", " vs ", results$comparison)


# ----- Plot

results <- results[results$padj < 0.05,]
results <- results[results$estimate != 0,]
results$l2fc <- log2(results$estimate)
results$log2upper <- log2(results$upper)
results$log2lower <- log2(results$lower)
results$point_size <- abs(results$l2fc) 

#custom label
custom_labels <- c(
  "BElike vs Naivelike" = "BElike vs **Naivelike**",
  "BElike vs Negative" = "BElike vs **Negative**",
  "Naivelike vs Negative" = "Naivelike vs **Negative**")


pdf("OCCAMS_EnrichmentDepletionMutSigComparison.pdf", w = 7, h = 6)
ggplot(results, aes(x = l2fc, y = Sig, color = -log10(padj), size = point_size)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#CC3311") +
  scale_color_gradient(low = "lightblue", high = "darkred", name = "log(q-value)") +
  scale_size_continuous(range = c(1, 10), name = "log2(OR)") +
  theme_minimal() +
  xlab("Log2 Odds Ratio") +
  ylab("Mutational Signature") +
  theme(axis.title.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "gray90", color = "black"), strip.text = ggtext::element_markdown()) + xlab("Log(Odds Ratio)") +
  facet_wrap(~comparison, nrow = 1, labeller = as_labeller(custom_labels)) +
  scale_x_continuous(breaks = c(-4,-2,0,2,4)) + expand_limits(x = c(-5, 5))
dev.off()

