setwd("/Users/jao/Desktop/MSc_Project/MMR")
library(dplyr)
library(wesanderson)
library(ggpubr)

# Load single base subsititution and their corresponding gene data
load("vep.snv_shorten.Rdata")
vep.snvs$Sample <- sub("_vs_.*", "", vep.snvs$Sample)

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))

# Select only primary samples
vep.snvs <- vep.snvs[vep.snvs$Sample %in% occams_results_df$Sample,]

# Extract only mmr genes
mmrgenes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "MSH3", "MLH3")
vep.snvs <- vep.snvs[vep.snvs$Gene %in% mmrgenes,]

# Count mutations per gene per sample, include zeros because they are biologically relevant

mmrMutCounts <- vep.snvs %>%
  group_by(Sample, Gene) %>%
  summarize(mutation_count = n())

mmrMutCounts_merged <- merge(mmrMutCounts, occams_results_df, by = "Sample")
mmrMutCounts_merged$Phenotype_Assigned <- factor(mmrMutCounts_merged$Phenotype_Assigned, levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                                      "TreatedLike.Sig17+", "Sig17-"))

my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

pdf("MMR_MutationsComparison.pdf", h = 6, w = 5)
ggboxplot(mmrMutCounts_merged, x = "Phenotype_Assigned", y = "mutation_count",
                  fill = "Phenotype_Assigned", palette = my_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) + ylab("Mutational Counts") +
  geom_pwc(method = "wilcox.test", label = "p.signif",
           p.adjust.method = "BH")
dev.off()


mmrMutCounts <- expand.grid(Sample = unique(occams_results_df$Sample), Gene = mmrgenes) %>%
  left_join(occams_results_df, by = "Sample") %>%
  left_join(mmrMutCounts, by = c("Sample", "Gene"))

mmrMutCounts$mutation_count[is.na(mmrMutCounts$mutation_count)] <- 0

mmrMutCounts$logcounts <- log1p(mmrMutCounts$mutation_count)

# ----- Plot 
my_colors <- c(
  "Sig17-" = wes_palette("GrandBudapest1")[2],
  "Sig17+" = "#CAB2D6",
  "BarrettsLike.Sig17+" = wes_palette("GrandBudapest1")[4])

pdf("MMR_Mutations3GroupComparison.pdf", h = 15, w = 10)
plot <- ggboxplot(mmrMutCounts, x = "Phenotype_Assigned", y = "logcounts",
         fill = "Phenotype_Assigned", palette = my_colors) +
  facet_wrap(~Gene, scales = "free",nrow = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) + ylab("Mutational Counts") +
  geom_pwc(method = "dunn_test", label = "p.adj.format",
           p.adjust.method = "BH",
           bracket.nudge.y = 0.2,
           step.increase = 0.15) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggadjust_pvalue(plot, p.adjust.method = "BH",
                label = "{p.adj.signif}", hide.ns = F)
dev.off()

# ---- Per sample

setwd("/Users/jao/Desktop/MSc_Project/MMR")
library(dplyr)
library(wesanderson)
library(ggpubr)

# Load single base subsititution and their corresponding gene data
load("vep.snv_shorten.Rdata")
vep.snvs$Sample <- sub("_vs_.*", "", vep.snvs$Sample)

# Load Sig17 Status
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))

# Select only primary samples
vep.snvs <- vep.snvs[vep.snvs$Sample %in% occams_results_df$Sample,]

# Extract only mmr genes
mmrgenes <- c("MLH1", "MSH2", "MSH6", "PMS1", "PMS2", "MSH3", "MLH3")
vep.snvs <- vep.snvs[vep.snvs$Gene %in% mmrgenes,]

# Count mutations per gene per sample, include zeros because they are biologically relevant

mmrMutCounts <- vep.snvs %>%
  group_by(Sample) %>%
  summarize(mutation_count = n())

mmrMutCounts_merged <- merge(mmrMutCounts, occams_results_df, by = "Sample")
mmrMutCounts_merged$Phenotype_Assigned <- factor(mmrMutCounts_merged$Phenotype_Assigned, levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                                                                    "TreatedLike.Sig17+", "Sig17-"))

my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

pdf("MMR_MutationsComparison_combined.pdf", h = 6, w = 5)
ggboxplot(mmrMutCounts_merged, x = "Phenotype_Assigned", y = "mutation_count",
          fill = "Phenotype_Assigned", palette = my_colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.x = element_blank()) + ylab("Mutational Counts") +
  geom_pwc(method = "dunn_test",
           p.adjust.method = "BH")
dev.off()

