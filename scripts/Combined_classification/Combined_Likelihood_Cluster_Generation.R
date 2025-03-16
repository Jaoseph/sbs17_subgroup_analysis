setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier")

# Load libraries

library(tidyverse)
library(deconstructSigs) # for context ordering
library(ggplot2)


load("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Naive_clust_mclust_meanCont.Rdata")
naive_dist <- mut.dists_mean

load("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/ICGC_Treated_clust_mclust_meanCont.Rdata")
treated_dist <- mut.dists_mean

load("/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_clust_mclust_meanCont.Rdata")
barrett_dist <- mut.dists_mean

rm(mut.dists_mean)

# Average Sig17- across all three  
negative_sig <- rbind(naive = naive_dist[rownames(naive_dist) %in% "Sig17-", ], 
                      treated = treated_dist[rownames(treated_dist) %in% "Sig17-", ],
                      barrett =barrett_dist[rownames(barrett_dist) %in% "Sig17-", ])
negative_sig <- rbind(mean = apply(negative_sig, 2, mean))


mut.dists_mean <- rbind(
  NaiveLike = naive_dist[rownames(naive_dist) %in% "Sig17+", ], 
  TreatedLike = treated_dist[rownames(treated_dist) %in% "Sig17+", ],
  BarrettsLike = barrett_dist[rownames(barrett_dist) %in% "Sig17+", ],
  Negative = negative_sig[rownames(negative_sig) %in% "mean"]
)

save(mut.dists_mean, file = 'Combined_clust_mclust_meanCont.Rdata')


# Plot prior clusters
mut.dists_mean.plot <- cbind(mut.dists_mean[,colnames(signatures.cosmic)],
                             mut.dists_mean[,97:ncol(mut.dists_mean)])
mut.dists_mean.plot$Phenotype <- rownames(mut.dists_mean)
mut.dists_mean.plot <- mut.dists_mean.plot %>%
  pivot_longer(cols = -Phenotype, names_to = 'Context', values_to = 'Contribution')
mut.dists_mean.plot$Phenotype <- factor(mut.dists_mean.plot$Phenotype,
                                        levels = sort(rownames(mut.dists_mean)))

# Sort indel context order
signatures.id83 <- read.table('COSMIC_v3.3_ID_GRCh37.txt', h=T)
mut.dists_mean.plot$Context <- factor(mut.dists_mean.plot$Context,
                                      levels = c(colnames(signatures.cosmic), signatures.id83$Type))
mut.dists_mean.plot$Type <- ifelse(
  grepl(pattern = '>', mut.dists_mean.plot$Context),
  'SBS','indel')

g_meanPlot <- ggplot(mut.dists_mean.plot, aes(x = Context, y = Contribution, fill = Type)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = 'top') +
  facet_wrap(~Phenotype, scales = 'free')

g_meanPlot
ggsave(filename = 'Combined_LikelihoodDistributionsMeans.pdf',
       plot = g_meanPlot, height = 5, width = 10)

