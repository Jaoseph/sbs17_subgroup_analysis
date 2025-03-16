setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/PCA")

#--- PCA Combined

library(ggfortify)
library(wesanderson)

# Load DeconstructSigs Mutational Signature Contributions for 721 samples
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names; primaries_sigs_complete <- primaries_sigs_complete[,-1]

#sig_phenotype <- primaries_sigs_complete$Phenotype_Assigned

primaries_sigs <- primaries_sigs_complete[, sapply(primaries_sigs_complete, is.numeric)]

#PCA
pca <- prcomp(primaries_sigs, scale = FALSE)

pdf("OCCAMS_PCA.pdf")
autoplot(pca, data = primaries_sigs_complete, color = "Phenotype_Assigned") + theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_color_manual(
    values = c(
      "Sig17-" = wes_palette("AsteroidCity3")[2],
      "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
      "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
      "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])) +
  ylab("Principal Component 2 (20.62%)") + xlab("Principal Component 1 (33.68%)") +
  labs(colour = "Signature Phenotype")
dev.off()
