# ------ Cohort-level dndscv - Prep (stratification of mutational data)
#setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")

# Load sig17 status annotations - Server
load("OCCAMS_CombinedresultsSummary.Rdata")

# Load Mutational Data - Server
load("mutational_OCCAMS.Rdata")

# Extract only Negative Samples - Server
negative <- mutational_OCCAMS[mutational_OCCAMS$SampleID %in% occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "Sig17-"],] #320 Sig17- Samples
negative <- negative[!negative$chr %in% c("hs37d5", "Un", "MT"),] 
length(unique(negative$SampleID)) #check
save(negative, file = "negativegroup.Rdata")

# Extract only BELike Samples - Server
BELike <- mutational_OCCAMS[mutational_OCCAMS$SampleID %in% occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"],] #101 BELike Samples
BELike <- BELike[!BELike$chr %in% c("hs37d5", "Un", "MT"),] 
length(unique(BELike$SampleID)) #check
save(BELike, file = "BELikegroup.Rdata")

# Extract only NaiveLike Samples - Server
NaiveLike <- mutational_OCCAMS[mutational_OCCAMS$SampleID %in% occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"],] #237 NaiveLike Samples
NaiveLike <- NaiveLike[!NaiveLike$chr %in% c("hs37d5", "Un", "MT"),] 
length(unique(NaiveLike$SampleID)) #check
save(NaiveLike, file = "NaiveLikegroup.Rdata")

# Extract only TreatedLike Samples - Server
TreatedLike <- mutational_OCCAMS[mutational_OCCAMS$SampleID %in% occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"],] #63 TreatedLike Samples
TreatedLike <- TreatedLike[!TreatedLike$chr %in% c("hs37d5", "Un", "MT"),] 
length(unique(TreatedLike$SampleID)) #check
save(TreatedLike, file = "TreatedLikegroup.Rdata")


