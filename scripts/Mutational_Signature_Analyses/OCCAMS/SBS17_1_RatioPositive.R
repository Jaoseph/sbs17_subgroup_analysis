# SBS17/SBS1 ratio comparison between BE and primaries

setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs")

# ----- Barretts Esophagus
load("Barretts_deconstructSigs.Rdata")
barretts_sigs_complete <- barretts_sigs_complete[,c("SBS1","SBS17a", "SBS17b")]

barretts_ratio <- barretts_sigs_complete/barretts_sigs_complete[,1] ; barretts_ratio <- barretts_ratio[,-1]
#barretts_ratio$tissue <- "Barretts"

load("/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata")
positiveBE <- rownames(ann)[ann$Phenotype == "Sig17+"]
barretts_ratio$status <- ifelse(rownames(barretts_ratio) %in% positiveBE, "BarrettsSig17+", "Sig17-")

# Remove Sig17-
barretts_ratio <- barretts_ratio[barretts_ratio$status != "Sig17-",]


# ----- Primaries
load("Primaries_deconstructSigs.Rdata")
primaries_sigs_complete <- primaries_sigs_complete[,c("SBS1","SBS17a", "SBS17b")]

primaries_ratio <- primaries_sigs_complete/primaries_sigs_complete[,1] ; primaries_ratio <- primaries_ratio[,-1]
#primaries_ratio$tissue <- "PrimaryTumour"

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)
primaries_ratio$status <- ifelse(rownames(primaries_ratio) %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                                 ifelse(rownames(primaries_ratio) %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                        ifelse(rownames(primaries_ratio) %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))

# Remove Sig17-
primaries_ratio <- primaries_ratio[primaries_ratio$status != "Sig17-",]

#------ Plot 

library(reshape2)
library(ggpubr)
library(wesanderson)

sig17 <- rbind(barretts_ratio, primaries_ratio)

sig17 <- melt(sig17, id.vars = "status", variable.name = "Signature", value.name = "Ratio")
sig17$l2fc <- log2(sig17$Ratio)

my_comparisons <- list(
  c("NaiveLike.Sig17+", "BarrettsLike.Sig17+"),
  c("NaiveLike.Sig17+", "TreatedLike.Sig17+"),
  c("NaiveLike.Sig17+", "BarrettsSig17+"),
  c("BarrettsLike.Sig17+", "TreatedLike.Sig17+"),
  c("BarrettsLike.Sig17+", "BarrettsSig17+"),
  c("TreatedLike.Sig17+", "BarrettsSig17+")
)
sig17$status <- factor(sig17$status, levels = c("BarrettsSig17+", "BarrettsLike.Sig17+",
                                              "NaiveLike.Sig17+", "TreatedLike.Sig17+"))

my_colors <- c(
  "BarrettsSig17+" = wes_palette("IsleofDogs1")[1],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])


#max_value <- max(sig17$Ratio, na.rm = TRUE) * 1.7
max_value <- max(sig17$l2fc, na.rm = TRUE) * 1.5

line_data <- data.frame(Signature = c("SBS17a", "SBS17b"),
                        yintercept = c(0,0))

pdf("Ratio_Analysis/OCCAMS_SBS17_1RatioBarrettsvsPrimaries.pdf", w = 8, h = 4)
ggviolin(sig17, x = "status", y = "l2fc", fill = "status", 
         palette = my_colors, add = "boxplot",facet.by = "Signature",
         add.params = list(width = 0.12, color = "black", fill = "white")) + 
  geom_pwc(method = "dunn_test", label = "p.adj.signif", hide.ns = TRUE,
           p.adjust.method = "BH",
           bracket.nudge.y = 0.02,
           step.increase = 0.12) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_blank(), axis.title.x =  element_blank(),
        axis.ticks.x = element_blank(), legend.text = element_text(size = 9)) + ylab("SBS17/SBS1 Ratio") +
  geom_hline(data = line_data, aes(yintercept = yintercept), linetype = "longdash", color = "red", size = 0.5)
dev.off()




  
  
  
  
  
