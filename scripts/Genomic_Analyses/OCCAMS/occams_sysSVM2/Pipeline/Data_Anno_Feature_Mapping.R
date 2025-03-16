# ---- Collate mutational and CNV data annotations

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
source("R/annotation_functions.R")

molecular_data = read_tsv("example_data/molecular_features_100samples.tsv")

# Load Small somatic mutation annotations
load("ssms_annotated.RData")

# Load CNV segment annotations
load("cnvs_annotated.Rdata")

# Collate into unified format for sysSVM2 : 718 samples..?
molecular_data = make_totalTable(
  ssms_annotated, cnvs_annotated, 
  canonical_drivers = "example_data/canonical_drivers.rds"
)

save(molecular_data, file = "molecular_data.Rdata")

# Load compendium of 25 of system-level properties of the genes
systemsLevel_data = read_tsv("example_data/systemsLevel_features_allGenes.tsv")

# Join the two tables to create the sysSVM2 input file
sysSVM2_input = inner_join(molecular_data, systemsLevel_data, by = "entrez")

save(sysSVM2_input, file = "sysSVM2_input.Rdata")


