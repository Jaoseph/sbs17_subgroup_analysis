setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Corr/New")
library(dplyr)

# Load DeconstructSigs Mutational Signature Contributions for 721 samples
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(primaries_sigs_complete > 0.05) > 0.2
primaries_sigs_complete <- primaries_sigs_complete[,keep]

primaries_sigs_complete <- primaries_sigs_complete[,names(primaries_sigs_complete)[grepl("SBS", names(primaries_sigs_complete))]]
primaries_sigs_complete <- primaries_sigs_complete[, c("SBS17a","SBS17b", "SBS18")]


#load classification annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")

# Remove Sig17- 
primaries_sigs_complete <- primaries_sigs_complete[primaries_sigs_complete$Phenotype_Assigned != "Sig17-",]

# Get statistics of how many samples have 0 SBS18 contribution in each group

test <- primaries_sigs_complete[primaries_sigs_complete$SBS18 == 0,]
table(test$Phenotype_Assigned)
#BarrettsLike.Sig17+    NaiveLike.Sig17+  TreatedLike.Sig17+ 
#  63                  33                  17 

# Remove SBS18 0 contribution samples 

primaries_sigs_complete <- primaries_sigs_complete[!primaries_sigs_complete$Row.names %in% test$Row.names,]

# ----- Plot

library(wesanderson)
my_colors <- c(
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

library(ggpubr)

pdf(file = "CorrScatterPlotSBS17b.pdf", h = 5, w = 6)
ggscatter(primaries_sigs_complete, x = "SBS17b", y= "SBS18", color = "Phenotype_Assigned", 
          add = "reg.line", conf.int = TRUE, palette = my_colors, text=element_text(size=10)) + facet_wrap(~Phenotype_Assigned, scale ="free")+
  stat_cor(method = "pearson")
dev.off()


sbs17a <- primaries_sigs_complete[primaries_sigs_complete$SBS17a > 0,]

pdf(file = "CorrScatterPlotSBS17a.pdf", h = 5, w = 6)
ggscatter(sbs17a, x = "SBS17a", y= "SBS18", color = "Phenotype_Assigned", 
          add = "reg.line", conf.int = TRUE, palette = my_colors, text=element_text(size=10)) + facet_wrap(~Phenotype_Assigned, scale ="free")+
  stat_cor(method = "pearson")
dev.off()
