setwd("/Users/jao/Desktop/MSc_Project/TCGA_deconstructSigs")

library(reshape2)
library(wesanderson)
library(ggpubr)
library(dplyr)
library(tidyr)

# ---------- Summary Statistics

# Load TCGA primary mut sig contributions
load("TCGA_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
tcga_sigs_complete[tcga_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(tcga_sigs_complete > 0.05) > 0.2
tcga_sigs_complete <- tcga_sigs_complete[,keep]


# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga <- subset(ann_tcga, select = Phenotype_Assigned)
tcga_sigs_complete <- merge(ann_tcga, tcga_sigs_complete, by = "row.names")
rownames(tcga_sigs_complete) <- tcga_sigs_complete$Row.names ; tcga_sigs_complete <- tcga_sigs_complete[,-1]

# Remove BarrettLike Sample
tcga_sigs_complete <- tcga_sigs_complete[tcga_sigs_complete$Phenotype_Assigned != "BarrettsLike.Sig17+", ,drop = FALSE]

# Factor sig 17 status and and SBS28
tcga_sigs_complete$Phenotype_Assigned <- factor(tcga_sigs_complete$Phenotype_Assigned, 
                                                     levels=c("NaiveLike.Sig17+", "Sig17-"))

summary_stats <- tcga_sigs_complete %>%
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

write.csv(summary_stats, file = "tcga_MutationalSigContributionSummary.csv", row.names = F)


# ---------- SBS
rm(list = ls())

# Load TCGA primary mut sig contributions
load("TCGA_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
tcga_sigs_complete[tcga_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(tcga_sigs_complete > 0.05) > 0.2
tcga_sigs_complete <- tcga_sigs_complete[,keep]

tcga_sigs_complete <- tcga_sigs_complete[, grepl('SBS', colnames(tcga_sigs_complete))]

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga <- subset(ann_tcga, select = Phenotype_Assigned)
tcga_sigs_complete <- merge(ann_tcga, tcga_sigs_complete, by = "row.names")
rownames(tcga_sigs_complete) <- tcga_sigs_complete$Row.names ; tcga_sigs_complete <- tcga_sigs_complete[,-1]

# Remove BarrettLike Sample
tcga_sigs_complete <- tcga_sigs_complete[tcga_sigs_complete$Phenotype_Assigned != "BarrettsLike.Sig17+", ,drop = FALSE]

# Factor sig 17 status and and SBS28
tcga_sigs_complete$Phenotype_Assigned <- factor(tcga_sigs_complete$Phenotype_Assigned, 
                                                levels=c("NaiveLike.Sig17+", "Sig17-"))

# Melt
melted_sigs_complete <- melt(tcga_sigs_complete)

# Reorder
melted_sigs_complete$variable <- factor(melted_sigs_complete$variable, 
                                        levels = c("SBS1", "SBS5", "SBS17b"))

my_comparisons <- list(
  c("NaiveLike.Sig17+", "Sig17-"))


my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4])


# Plot

pdf("tcga_SBSContributionDifferences.pdf", h = 8, w = 8)
ggviolin(melted_sigs_complete, x = "Phenotype_Assigned", y = "value",
                 fill = "Phenotype_Assigned", palette = my_colors,
                 add = "boxplot", add.params = list(width = 0.1, color = "black", fill = "white"))+
  xlab("")+ ylab("% Contribution") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~variable, 
             scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.signif",
           p.adjust.method = "BH",
           bracket.nudge.y = 0.5,
           step.increase = 0.3) + scale_y_continuous(,expand = expansion(mult = c(0, 0.3)))
dev.off()

# ---------- INDEL
rm(list = ls())

# Load TCGA primary mut sig contributions
load("TCGA_deconstructSigs.Rdata")
# Exclude if overall contribution less than 0.05 (5%)
tcga_sigs_complete[tcga_sigs_complete < 0.05] <- 0

# Exclude if only contributes in a small fraction of the samples (less than 20%)
keep <- colMeans(tcga_sigs_complete > 0.05) > 0.2
tcga_sigs_complete <- tcga_sigs_complete[,keep]

tcga_sigs_complete <- tcga_sigs_complete[, grepl('ID', colnames(tcga_sigs_complete))]

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga <- subset(ann_tcga, select = Phenotype_Assigned)
tcga_sigs_complete <- merge(ann_tcga, tcga_sigs_complete, by = "row.names")
rownames(tcga_sigs_complete) <- tcga_sigs_complete$Row.names ; tcga_sigs_complete <- tcga_sigs_complete[,-1]

# Remove BarrettLike Sample
tcga_sigs_complete <- tcga_sigs_complete[tcga_sigs_complete$Phenotype_Assigned != "BarrettsLike.Sig17+", ,drop = FALSE]

# Factor sig 17 status and and SBS28
tcga_sigs_complete$Phenotype_Assigned <- factor(tcga_sigs_complete$Phenotype_Assigned, 
                                                levels=c("NaiveLike.Sig17+", "Sig17-"))

# Melt
melted_sigs_complete <- melt(tcga_sigs_complete)

# Reorder
melted_sigs_complete$variable <- factor(melted_sigs_complete$variable, 
                                        levels = c("ID4", "ID5"))

my_comparisons <- list(
  c("NaiveLike.Sig17+", "Sig17-"))


my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4])



pdf("tcga_IDContributionDifferences.pdf", h = 8, w = 8)
ggviolin(melted_sigs_complete, x = "Phenotype_Assigned", y = "value",
                 fill = "Phenotype_Assigned", palette = my_colors,
                 add = "boxplot", add.params = list(width = 0.1, color = "black", fill = "white"))+
  xlab("")+ ylab("% Contribution") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) +
  facet_wrap(~variable, 
             scales = "free", nrow = 1)+ 
  geom_pwc(method = "dunn_test", label = "p.signif",
           p.adjust.method = "BH",
           bracket.nudge.y = 0.7,
           step.increase = 0.3) + scale_y_continuous(expand = expansion(mult = c(0, 0.3)))
dev.off()

