# Look at SBS18 and Sig17 correlations in BE

setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Corr/New")
library(dplyr)

# 2. Load annotation (ann)
load('/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata')

# Load BE DeconstructSigs Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata")
sigs_complete[sigs_complete < 0.05] <- 0
sigs_complete <- sigs_complete[, c("SBS17a","SBS17b", "SBS18")]

# join annotations and sigcontributions
sigs_complete <- merge(ann, sigs_complete, by = "row.names")

# Remove Sig17-
sigs_complete <- sigs_complete[sigs_complete$Phenotype != "Sig17-",]

# Remove SBS18 0 contribution samples 
test <- sigs_complete[sigs_complete$SBS18 == 0,]
sigs_complete <- sigs_complete[!sigs_complete$Row.names %in% test$Row.names,]

# Plot 

library(wesanderson)
"BarrettsSig17+" = wes_palette("IsleofDogs1")[1]

library(ggpubr)

pdf(file = "BECorrScatterPlotSBS17b.pdf", h = 5, w = 6)
ggscatter(sigs_complete, x = "SBS17b", y= "SBS18", color = "Phenotype", 
          add = "reg.line", conf.int = TRUE, palette = wes_palette("IsleofDogs1")[1], text=element_text(size=10)) +
  stat_cor(method = "pearson")
dev.off()


sbs17a <- sigs_complete[sigs_complete$SBS17a > 0,]

pdf(file = "BECorrScatterPlotSBS17a.pdf", h = 5, w = 6)
ggscatter(sbs17a, x = "SBS17a", y= "SBS18", color = "Phenotype", 
          add = "reg.line", conf.int = TRUE, palette = wes_palette("IsleofDogs1")[1], text=element_text(size=10))+
  stat_cor(method = "pearson")
dev.off()
