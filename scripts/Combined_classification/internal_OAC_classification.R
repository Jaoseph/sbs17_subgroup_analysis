#####
### Apply Sig17 classifier to OAC WGS internal cohort -
#####

#setwd("/Users/jao/Desktop/MSc_Project/Pan-cancer_Analysis/Analysis")
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

### Load All OAC WGS cohort data and tally mutation contributions

load("primaries.vcf.RData")
load("vcfs.extra.OAC.RData")
eac_uk <- read.table('simple_somatic_mutation.open.ESAD-UK.tsv.gz',
                     sep = '\t', header = TRUE)
eac_wgs <- eac_uk[eac_uk$sequencing_strategy == 'WGS', ]
eac_wgs <- eac_wgs[!duplicated(eac_wgs$icgc_mutation_id), ]
# Generate SNV Tally input

# Transform SNV data
print(" Transform SNV data")

primaries.vcf$SampleID <- sub("_vs_.*", "", primaries.vcf$Sample) 
oac_snv.vcf <- data.frame(
  Tumor_Sample_Barcode = primaries.vcf$SampleID,
  Chromosome = primaries.vcf$`#CHROM`,
  Start_Position = primaries.vcf$POS,
  End_Position = primaries.vcf$POS,
  Variant_Type = "SNP",
  Reference_Allele = primaries.vcf$REF,
  Tumor_Seq_Allele2 = primaries.vcf$ALT
)
oac_snv.vcf$Chromosome <- sub("chr","", oac_snv.vcf$Chromosome)

oac_extra.vcf <- data.frame(
  Tumor_Sample_Barcode = vcfs.extra.OAC$TumourID,
  Chromosome = vcfs.extra.OAC$Chromosome,
  Start_Position = vcfs.extra.OAC$Position,
  End_Position = vcfs.extra.OAC$Position,
  Variant_Type = "SNP",
  Reference_Allele = vcfs.extra.OAC$Ref,
  Tumor_Seq_Allele2 = vcfs.extra.OAC$Alt
)

icgc_snv.vcf <- data.frame(
  Tumor_Sample_Barcode = eac_wgs$submitted_sample_id,
  Chromosome = eac_wgs$chromosome,
  Start_Position = eac_wgs$chromosome_start,
  End_Position = eac_wgs$chromosome_end,
  Variant_Type = "SNP",
  Reference_Allele = eac_wgs$mutated_from_allele,
  Tumor_Seq_Allele2 = eac_wgs$mutated_to_allele
)

# Combine all WGS data
total_snv.vcf <- rbind(oac_snv.vcf,oac_extra.vcf,icgc_snv.vcf)
# Remove Duplicate Rows 
total_snv.vcf <- total_snv.vcf[!duplicated(total_snv.vcf),]

# Transform Indel Data
print(" Transform Indel data")

oac_indels <- read.delim("indels.OAC.txt")

# Split column into Chromosome, Position, Ref/Alt
chr_pos_change <- strsplit(as.character(oac_indels$V1), "_")
oac_indels$Position <- sapply(chr_pos_change, '[', 2)
oac_indels$Chromosome <- sapply(chr_pos_change, '[', 1)

# Split Ref/Alt 
indel_change <- sapply(chr_pos_change, '[', 3)
indel_split <- strsplit(as.character(indel_change), "/")
oac_indels$REF <- sapply(indel_split, '[', 1)
oac_indels$ALT <- sapply(indel_split, '[', 2)

oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample) 

oac_ID.vcf <- data.frame(
  Tumor_Sample_Barcode = oac_indels$SampleID,
  Chromosome = oac_indels$Chromosome,
  Start_Position = oac_indels$Position,
  End_Position = oac_indels$Position,
  Variant_Type = sapply(oac_indels$REF,
                        function(x) ifelse(x == "-", "INS", "DEL")),
  Reference_Allele = oac_indels$REF,
  Tumor_Seq_Allele2 = oac_indels$ALT
)

# Join SNV and Indel Mutation Data
vcf_oac_input <- rbind(total_snv.vcf, oac_ID.vcf)

# Transform into MAF format
mut.maf <- read_maf_minimal(vcf_oac_input)

# Generate tally mutation contributions (input data)

mt_tally.oac <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  use_syn = TRUE
)

oac_tally <- as.data.frame(cbind(mt_tally.oac$SBS_96, mt_tally.oac$ID_83))

# Extract Samples which have both SBS and ID mutational tallies
Tally_sum <- rowSums(oac_tally[,1:96])
oac_tally <- oac_tally[Tally_sum != 0,]

print("Saving oac_tally")
save(oac_tally, 
     file = 'oac_tally.Rdata')

### Run In-house Classifier on OAC WGS internal cohort

# Load tally mutation contributions (input data)
load('oac_tally.Rdata')

# Load Prior Pan-cancer cluster mean distributions (cluster_distributions = mut.dists_mean)
load('ICGC_clust20_mclust_meanCont.Rdata')

# Signature Phenotype Assignment (cluster_assign = pheno_assigned)
load('ICGC_IDnormalised_PhenotypeAnnotation_clust.Rdata')
pheno_assigned <- ann$Phenotype


# Likelihood function that aligns a dataset with the designated mean distributions
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
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
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  marginal.probs <- marginal.probs[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(marginal.probs)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x))))
  
  return(final_probs)
  
}

# Apply log-likelihood approach
results.OCCAMS_loglik <- likelihood_calc(input_data = oac_tally, 
                                       cluster_distributions = mut.dists_mean,
                                       cluster_assign = pheno_assigned)

colnames(results.OCCAMS_loglik) <- rownames(mut.dists_mean)

results.OCCAMS_df <- data.frame(
  Sample = rownames(results.OCCAMS_loglik),
  Phenotype_Assigned = apply(results.OCCAMS_loglik, 1,
                             function(x) names(x)[x==max(x)])
)
rownames(results.OCCAMS_df) <- NULL

save(results.OCCAMS_df, file = 'internal_OAC_resultsSummary.Rdata')


results.OCCAMS_df <- data.frame(
  Patient = rownames(results.OCCAMS_loglik),
  Phenotype_Assigned = apply(results.OCCAMS_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.OCCAMS_loglik, 1, max)
)
rownames(results.OCCAMS_df) <- NULL
save(results.OCCAMS_df, file = 'OCCAMS_resultsSummarywithPval.Rdata')

# Plot Sig17+ vs Sig17- Samples for Internal OAC WGS cohort

library(ggplot2)
library(wesanderson)
library(dplyr)

Plot_df <-results.OCCAMS_df %>%
  group_by(Phenotype_Assigned) %>%
  summarize(count = n())

Plot_df$Phenotype_Assigned <- factor(Plot_df$Phenotype_Assigned, levels = c("Sig17+", "Sig17-"))

Plot <- ggplot(Plot_df, aes(x = Phenotype_Assigned, y = count, fill = Phenotype_Assigned)) +
  geom_bar(stat = "identity") + 
  xlab("Mutational signature phenotype profile") +
  ylab("Number of Samples") +
  scale_fill_manual(values = wes_palette("Royal2", n = 2, type = "discrete")) +
  theme_minimal()

ggsave("Internal_OAC_Classification.pdf", plot = Plot, device = "pdf")

