setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/PCA")
library(ggfortify)

#--- SBS only PCA

# Load DeconstructSigs Mutational Signature Contributions for 721 samples
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names; primaries_sigs_complete <- primaries_sigs_complete[,-1]
primaries_sigs_complete <- primaries_sigs_complete[,names(primaries_sigs_complete) %in% c("Phenotype_Assigned",grep('SBS', colnames(primaries_sigs_complete), value = TRUE))]

primaries_sigs <- primaries_sigs_complete[, sapply(primaries_sigs_complete, is.numeric)]
primaries_sigs <- primaries_sigs[, grepl('SBS', colnames(primaries_sigs))]

# PCA
pca <- prcomp(primaries_sigs, scale = FALSE)

pdf("OCCAMS_SBS_PCA.pdf")
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
  ylab("Principal Component 2 (27.07%)") + xlab("Principal Component 1 (50.96%)") +
  labs(colour = "Signature Phenotype")
dev.off()

# Loadings

loadings <- as.data.frame(pca$rotation)

PC1 <- loadings %>%
  select(PC1) %>%
  arrange(desc(abs(PC1)))
PC1$Signatures <- rownames(PC1)
  
PC2 <- loadings %>%
  select(PC2) %>%
  arrange(desc(abs(PC2)))
PC2$Signatures <- rownames(PC2)

pdf("OCCAMS_SBS_PC1_Loading.pdf")
ggbarplot(PC1, x = "Signatures", y = "PC1", fill = "Signatures", legend = "none") +
  rotate_x_text(45) + theme(axis.text.x = element_text(size = 10), 
                            strip.background = element_rect(color = "black", fill = "white")) + xlab("Signatures") +
  ylab("Principal Component 1")
dev.off()

pdf("OCCAMS_SBS_PC2_Loading.pdf")
ggbarplot(PC2, x = "Signatures", y = "PC2", fill = "Signatures", legend = "none") +
  rotate_x_text(45) + theme(axis.text.x = element_text(size = 10), 
                            strip.background = element_rect(color = "black", fill = "white")) + xlab("Signatures") +
  ylab("Principal Component 2")
dev.off()

#--- SBS only PCA

# Load DeconstructSigs Mutational Signature Contributions for 721 samples
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names; primaries_sigs_complete <- primaries_sigs_complete[,-1]
primaries_sigs_complete <- primaries_sigs_complete[,names(primaries_sigs_complete) %in% c("Phenotype_Assigned",grep('SBS', colnames(primaries_sigs_complete), value = TRUE))]
primaries_sigs_complete <- primaries_sigs_complete[,!names(primaries_sigs_complete) %in% c("SBS17a","SBS17b")]

primaries_sigs <- primaries_sigs_complete[, sapply(primaries_sigs_complete, is.numeric)]
primaries_sigs <- primaries_sigs[, grepl('SBS', colnames(primaries_sigs))]

# PCA without SBS17a and SBS17b
pca <- prcomp(primaries_sigs, scale = FALSE)

pdf("OCCAMS_SBS_PCA_NoSBS17.pdf")
autoplot(pca, data = primaries_sigs_complete, color = "Phenotype_Assigned") + theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_color_manual(
    values = c(
      "Sig17-" =  wes_palette("AsteroidCity3")[2],
      "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
      "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
      "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])) +
  ylab("Principal Component 2 (27.07%)") + xlab("Principal Component 1 (50.96%)") +
  labs(colour = "Signature Phenotype")
dev.off()

# Loadings

loadings <- as.data.frame(pca$rotation)

PC1 <- loadings %>%
  select(PC1) %>%
  arrange(desc(abs(PC1)))
PC1$Signatures <- rownames(PC1)

PC2 <- loadings %>%
  select(PC2) %>%
  arrange(desc(abs(PC2)))
PC2$Signatures <- rownames(PC2)

pdf("OCCAMS_SBS_PC1_Loading.pdf")
ggbarplot(PC1, x = "Signatures", y = "PC1", fill = "Signatures", legend = "none") +
  rotate_x_text(45) + theme(axis.text.x = element_text(size = 10), 
                            strip.background = element_rect(color = "black", fill = "white")) + xlab("Signatures") +
  ylab("Principal Component 1")
dev.off()

pdf("OCCAMS_SBS_PC2_Loading.pdf")
ggbarplot(PC2, x = "Signatures", y = "PC2", fill = "Signatures", legend = "none") +
  rotate_x_text(45) + theme(axis.text.x = element_text(size = 10), 
                            strip.background = element_rect(color = "black", fill = "white")) + xlab("Signatures") +
  ylab("Principal Component 2")
dev.off()

