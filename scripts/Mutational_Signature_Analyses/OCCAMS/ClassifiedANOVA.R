setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/anova")

library(reshape2)
library(wesanderson)
library(ggpubr)
library(dplyr)
library(tidyr)

# ---------- Summary Statistics

# Load OCCAMS primary mut sig contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(primaries_sigs_complete > 0.05) > 0.2
primaries_sigs_complete <- primaries_sigs_complete[,keep]


# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample ; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)
primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names ; primaries_sigs_complete <- primaries_sigs_complete[,-1]

# Factor sig 17 status and and SBS28
primaries_sigs_complete$Phenotype_Assigned <- factor(primaries_sigs_complete$Phenotype_Assigned, 
                                                     levels=c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                              "TreatedLike.Sig17+", "Sig17-"))

summary_stats <- primaries_sigs_complete %>%
  group_by(Phenotype_Assigned) %>%
  summarise(
    across(everything(), list(Mean = ~mean(.x, na.rm = TRUE), 
                              Median = ~median(.x, na.rm = TRUE))))
summary_stats <- melt(summary_stats, id.vars = "Phenotype_Assigned")

summary_stats <- summary_stats %>%
  mutate(
    Type = ifelse(grepl("^SBS", summary_stats$variable), "SBS", "ID"),
    Number = as.numeric(gsub("[^0-9]", "", summary_stats$variable))
  ) %>%
  arrange(Type, Number) %>%
  select(-Type, -Number)

# Save
library(openxlsx)
write.xlsx(summary_stats, "MutationalSigContributionSummary.xlsx", sheetName = "Table S", rowNames=TRUE)


# ---------- SBS
rm(list = ls())

# Load OCCAMS primary mut sig contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(primaries_sigs_complete > 0.05) > 0.2
primaries_sigs_complete <- primaries_sigs_complete[,keep]

primaries_sigs_complete <- primaries_sigs_complete[, grepl('SBS', colnames(primaries_sigs_complete))]

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample ; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)
primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names ; primaries_sigs_complete <- primaries_sigs_complete[,-1]


# Factor sig 17 status
primaries_sigs_complete$Phenotype_Assigned <- factor(primaries_sigs_complete$Phenotype_Assigned, 
                                                     levels=c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                              "TreatedLike.Sig17+", "Sig17-"))

# Melt
melted_sigs_complete <- melt(primaries_sigs_complete)

# Reorder
melted_sigs_complete$variable <- factor(melted_sigs_complete$variable, 
                                        levels = c("SBS1", "SBS5", "SBS17a", "SBS17b", "SBS18", "SBS40"))

melted_sigs_complete$variable <- factor(melted_sigs_complete$variable, 
                                        levels = c("SBS1", "SBS5", "SBS17a", "SBS17b", "SBS18", "SBS40", "ID1", "ID2", "ID4", "ID5", "ID9", "ID14"))

my_comparisons <- list(
  c("NaiveLike.Sig17+", "BarrettsLike.Sig17+"),
  c("NaiveLike.Sig17+", "TreatedLike.Sig17+"),
  c("NaiveLike.Sig17+", "Sig17-"),
  c("BarrettsLike.Sig17+", "TreatedLike.Sig17+"),
  c("BarrettsLike.Sig17+", "Sig17-"),
  c("TreatedLike.Sig17+", "Sig17-")
)


my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])



# help visualization
melted_sigs_complete <- melted_sigs_complete %>%
  mutate(value = ifelse(variable == "SBS17a" & value == 0 & Phenotype_Assigned == "Sig17-", 0.000001, value))


# Plot
pdf("SBSContributionDifferences.pdf", h = 10, w = 6)
plot <- ggviolin(melted_sigs_complete, x = "Phenotype_Assigned", y = "value",
                 fill = "Phenotype_Assigned", palette = my_colors,
                 add = "boxplot", add.params = list(width = 0.1, color = "black", fill = "white"))+
  xlab("")+ ylab("Mutational Signature Contribution") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~variable, 
             scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.adj.format",
           p.adjust.method = "BH",
           bracket.nudge.y = 0.1,
           step.increase = 0.3) + scale_y_continuous(,expand = expansion(mult = c(0, 0.1)))
ggadjust_pvalue(plot, p.adjust.method = "BH",
                label = "{p.adj.signif}", hide.ns = T)
dev.off()

# ---------- INDEL
rm(list = ls())

# Load OCCAMS primary mut sig contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(primaries_sigs_complete > 0.05) > 0.2
primaries_sigs_complete <- primaries_sigs_complete[,keep]

primaries_sigs_complete <- primaries_sigs_complete[, grepl('ID', colnames(primaries_sigs_complete))]

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
rownames(occams_results_df) <- occams_results_df$Sample ; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)
primaries_sigs_complete <- merge(occams_results_df, primaries_sigs_complete, by = "row.names")
rownames(primaries_sigs_complete) <- primaries_sigs_complete$Row.names ; primaries_sigs_complete <- primaries_sigs_complete[,-1]


# Factor sig 17 status 
primaries_sigs_complete$Phenotype_Assigned <- factor(primaries_sigs_complete$Phenotype_Assigned, 
                                                     levels=c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                              "TreatedLike.Sig17+", "Sig17-"))

# Melt
melted_sigs_complete <- melt(primaries_sigs_complete)

# Reorder
melted_sigs_complete$variable <- factor(melted_sigs_complete$variable, 
                                        levels = c("ID1", "ID2", "ID4", "ID5", "ID9", "ID14"))

my_comparisons <- list(
  c("NaiveLike.Sig17+", "BarrettsLike.Sig17+"),
  c("NaiveLike.Sig17+", "TreatedLike.Sig17+"),
  c("NaiveLike.Sig17+", "Sig17-"),
  c("BarrettsLike.Sig17+", "TreatedLike.Sig17+"),
  c("BarrettsLike.Sig17+", "Sig17-"),
  c("TreatedLike.Sig17+", "Sig17-")
)

my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

# Remove NS
melted_sigs_complete <- melted_sigs_complete %>%
  filter(!variable %in% c("ID10", "ID5"))


pdf("IDContributionDifferences.pdf", h = 10, w = 6)
plot <- ggviolin(melted_sigs_complete, x = "Phenotype_Assigned", y = "value",
                 fill = "Phenotype_Assigned", palette = my_colors,
                 add = "boxplot", add.params = list(width = 0.1, color = "black", fill = "white"))+
  xlab("")+ ylab("Mutational Signature Contribution") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~variable, 
             scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.adj.format",
           p.adjust.method = "BH",
           bracket.nudge.y = 0.1,
           step.increase = 0.3) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggadjust_pvalue(plot, p.adjust.method = "BH",
                label = "{p.adj.signif}", hide.ns = T)
dev.off()

