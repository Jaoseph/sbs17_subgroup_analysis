setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs")

load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")

primaries_sigs_complete <- primaries_sigs_complete[, c("SBS17a", "SBS17b")]
primaries_sigs_complete$Sample <- rownames(primaries_sigs_complete)


library(ggpubr)
library(wesanderson)


pdf("SBS17A_SBS17B_cor.pdf")
ggscatter(primaries_sigs_complete, x = "SBS17a", y = "SBS17b", color = "grey",
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, 
          cor.method = "pearson", xlab = "SBS17a Contribution", 
          ylab = "SBS17b Contribution", add.params = list(color = "#FC4E07"))
dev.off()


samples <- read.csv("/Users/jao/Desktop/MSc_Project/classifier_stratification/classifier_stratified_samples.csv")

primaries_sigs_complete <- merge(primaries_sigs_complete, samples, by = "Sample")

pdf("SBS17A_SBS17B_cor_color.pdf")
ggscatter(primaries_sigs_complete, x = "SBS17a", y = "SBS17b", color = "assignment",
          add = "reg.line", conf.int = TRUE, cor.coef = TRUE, 
          cor.method = "pearson", xlab = "SBS17a Contribution", 
          ylab = "SBS17b Contribution", add.params = list(color = "#FC4E07")) +
  scale_color_manual(values = c("BELike.Sig17+" = wes_palette("GrandBudapest1")[1],
                                "Sig17-" = wes_palette("GrandBudapest1")[2],
                                "TreatedLike.Sig17+" = wes_palette("GrandBudapest1")[3],
                                "NaiveLike.Sig17+" = wes_palette("GrandBudapest1")[4]))
dev.off()
