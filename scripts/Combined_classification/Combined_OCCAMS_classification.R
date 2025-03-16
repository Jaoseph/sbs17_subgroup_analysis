##### Apply Combined classifier to OCCAMS samples 
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier")

# ----- Loading.. 

# Load tally mutation contributions (input data)
load('oac_tally.Rdata')

# Load samplesID data, extract only primaries with data
#load("annotation.sampleIDs.RData")
#samples <- annotation.sampleIDs$TumourID[annotation.sampleIDs$Category == "PrimaryTumour"]
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
samples <- rownames(primaries_sigs_complete)

oac_tally <- oac_tally[rownames(oac_tally) %in% samples,]

# Load combined cluster mean distributions (cluster_distributions = mut.dists_mean)
load('/Users/jao/Desktop/MSc_Project/Combined_Classifier/Combined_clust_mclust_meanCont.Rdata')
rownames(mut.dists_mean) <- c("NaiveLike.Sig17+", "TreatedLike.Sig17+", "BarrettsLike.Sig17+", "Sig17-")

# ----- Attempt 2: EAC-specific Priors

# Barrets
load('/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata')
BEpheno_assigned <- ann$Phenotype
BEPositive <- sum(BEpheno_assigned == "Sig17+")
#BEmarginal.probs <- table(BEpheno_assigned)/length(BEpheno_assigned)

# Naive 
load('/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Naive_IDnormalised_PhenotypeAnnotation_clust.Rdata')
ann <- ann[ann$project_code == "ESAD",]
Naivepheno_assigned <- ann$Phenotype
NaivePositive <- sum(Naivepheno_assigned == "Sig17+")
#Naivemarginal.probs <- table(Naivepheno_assigned)/length(Naivepheno_assigned)

# Treated 
load('/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/ICGC_Treated_IDnormalised_PhenotypeAnnotation_clust.Rdata')
ann <- ann[ann$project_code == "ESAD",]
Treatedpheno_assigned <- ann$Phenotype
TreatedPositive <- sum(Treatedpheno_assigned == "Sig17+")
#Treatedmarginal.probs <- table(Treatedpheno_assigned)/length(Treatedpheno_assigned)

# Total Number of EAC samples (including BE) across all three
total_samples <- length(Treatedpheno_assigned) + length(Naivepheno_assigned) + length(BEpheno_assigned)

global_prior <- c(
  "NaiveLike_Sig17+" = NaivePositive/total_samples,
  "TreatedLike_Sig17+" = TreatedPositive/total_samples,
  "BarrettsLike_Sig17+" = BEPositive/total_samples,
  "Sig17-" = (total_samples-(NaivePositive+TreatedPositive+BEPositive))/total_samples
)

names(global_prior) <- c("NaiveLike.Sig17+", "TreatedLike.Sig17+", "BarrettsLike.Sig17+", "Sig17-")


# ----- Likelihood Function Load
#(cluster_assign = pheno_assigned)

# Likelihood function that aligns a dataset with the designated mean distributions
likelihood_calc <- function(input_data, cluster_distributions, prior) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   For now, this function does not limit mutation types: SBS and indels will be included
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  print('Calculating log-likelihoods...')
  for (i in 1:nrow(input_data)) {
    print(paste0('Calculating log likelihoods for sample ', i, ' of ', nrow(input_data), ': ', rownames(input_data)[i]))
    log_likelihoods[i, ] <- apply(cluster_distributions, 1, 
                                  function(x) sum(log10(x) * input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  prior <- prior
  prior <- prior[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(prior)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x))))
  
  return(final_probs)
  
}

# ----- Apply log-likelihood approach

results_loglik <- likelihood_calc(input_data = oac_tally, 
                                  cluster_distributions = mut.dists_mean,
                                  prior = global_prior)

colnames(results_loglik) <- rownames(mut.dists_mean)

occams_results_df <- data.frame(
  Sample = rownames(results_loglik),
  Phenotype_Assigned = apply(results_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results_loglik, 1, max)
)

rownames(occams_results_df) <- NULL
save(occams_results_df, file = 'OCCAMS_CombinedresultsSummary.Rdata')


# Save
library(openxlsx)
write.xlsx(occams_results_df, "OCCAMS_Combinedclassification_annotation.xlsx", sheetName = "Table S", rowNames=TRUE)

# ----- Plot

library(wesanderson)
library(dplyr)
library(ggpubr)

Plot_df <- occams_results_df %>%
  group_by(Phenotype_Assigned) %>%
  summarize(count = n()) %>%
  mutate(proportion = count / sum(count))

Plot_df$Phenotype_Assigned <- factor(Plot_df$Phenotype_Assigned, levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+",
                                                                            "TreatedLike.Sig17+", "Sig17-"))

Plot <- ggplot(Plot_df, aes(x = Phenotype_Assigned, y = count, fill = Phenotype_Assigned)) +
  geom_bar(stat = "identity") + 
  xlab("Mutational signature phenotype profile") +
  ylab("Number of Samples") +
  scale_fill_manual(values = c(
    "Sig17-" =  wes_palette("AsteroidCity3")[2],
    "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
    "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
    "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          legend.position="none")

Plot

ggsave("OCCAMS_CombinedClassification.pdf", plot = Plot, device = "pdf")
