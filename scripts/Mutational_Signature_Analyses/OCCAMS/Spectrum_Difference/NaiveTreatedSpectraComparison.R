##### Comparing mutational signature contribution patterns

setwd("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison")
library(dplyr)
library(ggpubr)

# ----- Naive vs Treated 

## Load Mean Representative Mutational Spectra

#Naive
load("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Naive_clust_mclust_meanCont.Rdata")
naive.mut.dists_mean <- mut.dists_mean
naive.mut.dists_mean <- naive.mut.dists_mean[!rownames(naive.mut.dists_mean) %in% "Sig17-",]

#Treated
load("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/ICGC_Treated_clust_mclust_meanCont.Rdata")
treated.mut.dists_mean <- mut.dists_mean
treated.mut.dists_mean <- treated.mut.dists_mean[!rownames(treated.mut.dists_mean) %in% "Sig17-",]

rm(mut.dists_mean)

## Calculate Mean Representative Mutational Spectra Difference

NT_diff <- naive.mut.dists_mean - treated.mut.dists_mean

NF_diff_plot <- cbind(NT_diff[,colnames(signatures.cosmic)], NT_diff[,97:ncol(NT_diff)])
NF_diff_plot$Phenotype <- rownames(NT_diff)
NF_diff_plot <- NF_diff_plot %>%
  filter(Phenotype == "Sig17+") %>%
  pivot_longer(cols = -Phenotype, names_to = 'Context', values_to = 'Contribution') %>%
  select(Context, Contribution)

signatures.id83 <- read.table('/Users/jao/Desktop/MSc_Project/Barretts_Classifier/MutationalSignatureReferences/COSMIC_v3.3_ID_GRCh37.txt', h=T)
NF_diff_plot$Context <- factor(NF_diff_plot$Context,
                               levels = c(colnames(signatures.cosmic), signatures.id83$Type))
NF_diff_plot$Type <- ifelse(grepl(pattern = '>', NF_diff_plot$Context),'SBS','indel')

# Remove non-significant mutational contexts
significance <- read.csv("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/Statistical_testing/NaiveTreated_signifiDiffSpectrum.csv")$Context
NF_diff_plot$Contribution <- ifelse(NF_diff_plot$Context %in% significance, NF_diff_plot$Contribution, 0)


## Plot SBS Mean Representative Mutational Spectra Difference

NF_diff_plot_SBS <- NF_diff_plot[NF_diff_plot$Type == 'SBS',]
NF_diff_plot_SBS$Substitution <- gsub(".*\\[(.*?)\\].*", "\\1", NF_diff_plot_SBS$Context)

NF_diff_plot_SBS$triplet <- NF_diff_plot_SBS$Context
NF_diff_plot_SBS$triplet <- gsub("\\[C>A\\]", "C", NF_diff_plot_SBS$triplet)
NF_diff_plot_SBS$triplet <- gsub("\\[C>G\\]", "C", NF_diff_plot_SBS$triplet)
NF_diff_plot_SBS$triplet <- gsub("\\[C>T\\]", "C", NF_diff_plot_SBS$triplet)
NF_diff_plot_SBS$triplet <- gsub("\\[T>A\\]", "T", NF_diff_plot_SBS$triplet)
NF_diff_plot_SBS$triplet <- gsub("\\[T>C\\]", "T", NF_diff_plot_SBS$triplet)
NF_diff_plot_SBS$triplet <- gsub("\\[T>G\\]", "T", NF_diff_plot_SBS$triplet)

sbs_meanPlot <- ggbarplot(NF_diff_plot_SBS, "Context", "Contribution", fill = "Substitution", 
                        color = "Substitution", palette = c("C>A" = "#0fbbe7","C>G" = "#0f0f0f", "C>T" = "#db2d28",
                                                            "T>A" = "#cdc9ca","T>C" = "#a3ce61","T>G" = "#ecc6c2", legend = "left")) + 
  rotate_x_text(90) + theme(axis.text.x = element_text(size = 8)) +
  ylab("Contribution Difference") + scale_x_discrete(labels = NF_diff_plot_SBS$triplet)

sbs_meanPlot
ggsave(filename = 'Naive-Treated_SBS.pdf',
       plot = sbs_meanPlot, height = 5, width = 10)

## Plot Indels Mean Representative Mutational Spectra Difference

NF_diff_plot_indel <- NF_diff_plot[NF_diff_plot$Type == 'indel',]
NF_diff_plot_indel$Context <- as.character(NF_diff_plot_indel$Context)

NF_diff_plot_indel$mutation_type <- ifelse(grepl(".*Del:C.*", NF_diff_plot_indel$Context) | grepl(".*Del:T.*", NF_diff_plot_indel$Context), "1bp Deletion",
                                           ifelse(grepl(".*Ins:C.*", NF_diff_plot_indel$Context) | grepl(".*Ins:T.*", NF_diff_plot_indel$Context), "1bp Insertion",
                                                  ifelse(grepl(".*Del:R.*", NF_diff_plot_indel$Context), ">1bp Deletion at Repeats (Deletion Length)",
                                                         ifelse(grepl(".*Ins:R.*", NF_diff_plot_indel$Context), ">1bp Insertion at Repeats (Insertion Length)", "Microhomology (Deletion Length)"))))
NF_diff_plot_indel$mutation_type <- factor(NF_diff_plot_indel$mutation_type, levels = c("1bp Deletion", "1bp Insertion", ">1bp Deletion at Repeats (Deletion Length)", ">1bp Insertion at Repeats (Insertion Length)",
                                                                                        "Microhomology (Deletion Length)"))


NF_diff_plot_indel$NContext <- ifelse(grepl(".*Del:C.*", NF_diff_plot_indel$Context), "C",
                                      ifelse(grepl(".*Del:T.*", NF_diff_plot_indel$Context), "T",
                                             ifelse(grepl(".*Ins:C.*", NF_diff_plot_indel$Context), "C",
                                                    ifelse(grepl(".*Ins:T.*", NF_diff_plot_indel$Context), "T",
                                                           ifelse(grepl("2:Del:R.*", NF_diff_plot_indel$Context), "2",
                                                                  ifelse(grepl("3:Del:R.*", NF_diff_plot_indel$Context), "3",
                                                                         ifelse(grepl("4:Del:R.*", NF_diff_plot_indel$Context), "4",
                                                                                ifelse(grepl("5:Del:R.*", NF_diff_plot_indel$Context), "5+",
                                                                                       ifelse(grepl("2:Ins:R.*", NF_diff_plot_indel$Context), "2",
                                                                                              ifelse(grepl("3:Ins:R.*", NF_diff_plot_indel$Context), "3",
                                                                                                     ifelse(grepl("4:Ins:R.*", NF_diff_plot_indel$Context), "4",
                                                                                                            ifelse(grepl("5:Ins:R.*", NF_diff_plot_indel$Context), "5+",
                                                                                                                   ifelse(grepl("2:Del:M.*", NF_diff_plot_indel$Context), "2",
                                                                                                                          ifelse(grepl("3:Del:M.*", NF_diff_plot_indel$Context), "3",
                                                                                                                                 ifelse(grepl("4:Del:M.*", NF_diff_plot_indel$Context), "4",
                                                                                                                                        "5+"))))))))))))))) 
NF_diff_plot_indel$Repeats <- sapply(NF_diff_plot_indel$Context, 
                                     function(x) substr(x, 9,9))
NF_diff_plot_indel$Repeats <- as.numeric(NF_diff_plot_indel$Repeats)
NF_diff_plot_indel <- NF_diff_plot_indel %>%
  mutate(Repeats = ifelse(row_number() <= 24, Repeats + 1, Repeats))
NF_diff_plot_indel$Repeats <- as.character(NF_diff_plot_indel$Repeats)
NF_diff_plot_indel$Context_Repeats <- paste(NF_diff_plot_indel$NContext, NF_diff_plot_indel$Repeats, sep = ":")

indel_meanPlot <- ggbarplot(NF_diff_plot_indel, x = "Context_Repeats", y = "Contribution", 
                            fill = "Context", color = "Context", legend = "none", 
                            facet.by = "mutation_type", nrow = 1, scales = "free_x", 
                            palette = ifelse(grepl(".*Del:C.*", NF_diff_plot_indel$Context), "#fdbe6e",
                                             ifelse(grepl(".*Del:T.*", NF_diff_plot_indel$Context), "#fe7e01",
                                                    ifelse(grepl(".*Ins:C.*", NF_diff_plot_indel$Context), "#acdc88",
                                                           ifelse(grepl(".*Ins:T.*", NF_diff_plot_indel$Context), "#38a02f",
                                                                  ifelse(grepl("2:Del:R.*", NF_diff_plot_indel$Context), "#facab4",
                                                                         ifelse(grepl("3:Del:R.*", NF_diff_plot_indel$Context), "#fe876f",
                                                                                ifelse(grepl("4:Del:R.*", NF_diff_plot_indel$Context), "#ed4633",
                                                                                       ifelse(grepl("5:Del:R.*", NF_diff_plot_indel$Context), "#bc191c",
                                                                                              ifelse(grepl("2:Ins:R.*", NF_diff_plot_indel$Context), "#d0e0ed",
                                                                                                     ifelse(grepl("3:Ins:R.*", NF_diff_plot_indel$Context), "#94c3e1",
                                                                                                            ifelse(grepl("4:Ins:R.*", NF_diff_plot_indel$Context), "#4c95c8",
                                                                                                                   ifelse(grepl("5:Ins:R.*", NF_diff_plot_indel$Context), "#1662ac",
                                                                                                                          ifelse(grepl("2:Del:M.*", NF_diff_plot_indel$Context), "#e2dff0",
                                                                                                                                 ifelse(grepl("3:Del:M.*", NF_diff_plot_indel$Context), "#b4b7d8",
                                                                                                                                        ifelse(grepl("4:Del:M.*", NF_diff_plot_indel$Context), "#8584bd",
                                                                                                                                               "#614099")))))))))))))))) + 
  rotate_x_text(90) + theme(axis.text.x = element_text(size = 8), 
                            strip.background = element_rect(color = "black", fill = "white")) + xlab("Context")

indel_meanPlot

ggsave(filename = 'Naive-Treated_Indels.pdf',
       plot = indel_meanPlot, height = 7, width = 16)
