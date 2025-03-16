#---- OCCAMS overall SNV, INDEL, CNV alteration comparison
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/genomic_changes")

library(dplyr)
library(wesanderson)
library(ggpubr)

# ----- Load alteration data

# Load CNV data (AMP and DEL)
oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$SampleID <- sub("_vs_.*", "", oac_cna$sample)

oac_del <- oac_cna[oac_cna$CNchange == "DEL", c("CNchange", "SampleID")]
oac_del$Alteration <- "Deletion"
names(oac_del) <- c("V7", "SampleID", "Alteration")
oac_amp <- oac_cna[oac_cna$CNchange == "AMP", c("CNchange", "SampleID")]
oac_amp$Alteration <- "Amplification"
names(oac_amp) <- c("V7", "SampleID", "Alteration")

# Load Nonsynonymous SNV data
oac_SNVs.nonsynonymous <- read.delim("SNVs.nonsynonymous.OAC.txt")
oac_SNVs.nonsynonymous$SampleID <- sub("_vs_.*", "", oac_SNVs.nonsynonymous$Sample)
oac_SNVs.nonsynonymous <- oac_SNVs.nonsynonymous[,c("V7", "SampleID")] 
oac_SNVs.nonsynonymous$Alteration <- "SNV"

# Load Indels data
oac_indels <- read.delim("indels.OAC.txt")
oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample)  
oac_indels <- oac_indels[,c("V7", "SampleID")] 
oac_indels$Alteration <- "Indel"


# Combine all data into one 
combined_oac <- rbind(oac_del, oac_amp, oac_indels, oac_SNVs.nonsynonymous)

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

# Extract only the 721 primary samples
combined_oac <- combined_oac[combined_oac$SampleID %in% rownames(occams_results_df),]

# ----- Count Alterations

# Amplification
amp_counts <- combined_oac %>%
  filter(Alteration == "Amplification") %>%
  group_by(SampleID) %>%
  summarise(count = n())
amp_counts$Group <- ifelse(amp_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                                 ifelse(amp_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                        ifelse(amp_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))
amp_counts <- amp_counts %>%
  mutate(log_count = log10(count))
amp_counts$Alteration <- "AMP"

# Deletion
del_counts <- combined_oac %>%
  filter(Alteration == "Deletion") %>%
  group_by(SampleID) %>%
  summarise(count = n())
del_counts$Group <- ifelse(del_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                           ifelse(del_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                  ifelse(del_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))
del_counts <- del_counts %>%
  mutate(log_count = log10(count))
del_counts$Alteration <- "DEL"

# Indels
indel_counts <- combined_oac %>%
  filter(Alteration == "Indel") %>%
  group_by(SampleID) %>%
  summarise(count = n())
indel_counts$Group <- ifelse(indel_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                           ifelse(indel_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                  ifelse(indel_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))
indel_counts <- indel_counts %>%
  mutate(log_count = log10(count))
indel_counts$Alteration <- "Indel"

# Non-synoymous mutations
snv_counts <- combined_oac %>%
  filter(Alteration == "SNV") %>%
  group_by(SampleID) %>%
  summarise(count = n())
snv_counts$Group <- ifelse(snv_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                             ifelse(snv_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                    ifelse(snv_counts$SampleID %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))
snv_counts <- snv_counts %>%
  mutate(log_count = log10(count))
snv_counts$Alteration <- "SNV"

# Collate all alteration counts
alterations_counts <- rbind(snv_counts, indel_counts, amp_counts, del_counts)
alterations_counts$Group <- factor(alterations_counts$Group, levels = c("BarrettsLike.Sig17+","NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-"))
alterations_counts$Alteration <- factor(alterations_counts$Alteration, levels = c("SNV", "Indel", "AMP", "DEL"))

# ----- Plot

my_comparisons <- list(
  c("NaiveLike.Sig17+", "BarrettsLike.Sig17+"),
  c("NaiveLike.Sig17+", "TreatedLike.Sig17+"),
  c("NaiveLike.Sig17+", "Sig17-"),
  c("BarrettsLike.Sig17+", "TreatedLike.Sig17+"),
  c("BarrettsLike.Sig17+", "Sig17-"),
  c("TreatedLike.Sig17+", "Sig17-"))

my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])


pdf("OCCAMS_overall_alterationComparison.pdf")
ggviolin(alterations_counts, x = "Group", y = "log_count", fill = "Group",
         palette = my_colors, add = "boxplot", add.params = list(width = 0.2, color = "black", fill = "white")) + 
  facet_wrap(alterations_counts$Alteration, scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.adj", hide.ns = FALSE,
           p.adjust.method = "BH",
           bracket.nudge.y = 0.1,
           step.increase = 0.1) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(), legend.text = element_text(size = 8)) +
  labs(colour = "Signature Phenotype") + ylab("Log10 Alteration Counts") + xlab("")
dev.off()

# Collate SNV and Indels only alteration counts
snv_indel_counts <- rbind(snv_counts, indel_counts)
snv_indel_counts$Group <- factor(snv_indel_counts$Group, levels = c("BarrettsLike.Sig17+","NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-"))
snv_indel_counts$Alteration <- factor(snv_indel_counts$Alteration, levels = c("SNV", "Indel"))

pdf("OCCAMS_SNVIndel_alterationComparison.pdf", w = 7, h = 12)
ggviolin(snv_indel_counts, x = "Group", y = "log_count", fill = "Group",
         palette = my_colors, add = "boxplot", add.params = list(width = 0.2, color = "black", fill = "white"),
         legend.title = "Signature Phenotype") + 
  facet_wrap(snv_indel_counts$Alteration, scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.adj.signif", hide.ns = FALSE,
           p.adjust.method = "BH",
           bracket.nudge.y = 0.2,
           step.increase = 0.1) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.text = element_text(size = 8), legend.title = element_text(size = 8), axis.text.y = element_text(size = 8)) +
  ylab("Log10 Alteration Counts")
dev.off()
