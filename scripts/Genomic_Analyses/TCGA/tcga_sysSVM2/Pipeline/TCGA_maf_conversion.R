# ----- Convert TCGA mutational and copy number data into MAF format to generate ANNOVAR input for sysSVM2 - DEPRECIATED

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA")

# Load libraries
library(deconstructSigs)
library(TCGAbiolinks)
library(sigminer)
library(maftools)
library(tidyverse)

# ---- Generate TCGA mutation data
# Load TCGA data 
query <- GDCquery(
  project = "TCGA-ESCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
#GDCdownload(query)
tcga_mutations <- GDCprepare(query = query)

# Exclude mutations from Esophageal Squamous Cell Carcinoma
tcga_annotation <- read.delim('esca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep="\t")
tcga_annotation <- tcga_annotation[tcga_annotation$Cancer.Type.Detailed == "Esophageal Adenocarcinoma",]
tcga_eac_ID_index <- tcga_annotation$Patient.ID
tcga_mutations$Patient_ID <- sapply(tcga_mutations$Tumor_Sample_Barcode,
                                    function(x) substr(x,1,12))
tcga_mutations <- tcga_mutations[which(tcga_mutations$Patient_ID %in% tcga_eac_ID_index),]

# Exclude mutations from non-primary tumours
tcga_mutations$sample_type_code <- sapply(tcga_mutations$Tumor_Sample_Barcode,
                                          function(x) substr(x,14,15))
tcga_mutations <- tcga_mutations[tcga_mutations$sample_type_code == '01', ]


# generate MAF format 
mut.maf <- read.maf(tcga_mutations, isTCGA = TRUE,
                    vc_nonSyn = names(table(tcga_mutations$Variant_Classification)))

write.mafSummary(maf = mut.maf, basename = 'esca')
