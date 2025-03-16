setwd("/Users/jao/Desktop/MSc_Project/TCGA_deconstructSigs")

# Load libraries
library(deconstructSigs)
library(TCGAbiolinks)
library(sigminer)
library(maftools)

# ---- Generate TCGA Mutation Tally
# Load TCGA data and tally mutation contributions
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

mut.maf <- read.maf(tcga_mutations, isTCGA = TRUE,
                    vc_nonSyn = names(table(tcga_mutations$Variant_Classification)))

# Tally 
mt_tally.tcga <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
  mode = 'ALL',
  use_syn = TRUE
)
tcga_tally <- as.data.frame(cbind(mt_tally.tcga$SBS_96, mt_tally.tcga$ID_83))


# ---- Format DeconstructSigs input

# Load selected EAC signature contributions to include
#   Only signatures appearing in >1% samples will be included
barretts_final.sigs <- read.table('EAC_sigProfilerCont.txt', header = TRUE)
barretts_final.sigs <- barretts_final.sigs[barretts_final.sigs$Prop_present > 0.01, ]
#sigs.sbs <- barretts_final.sigs$Sigs[grepl(pattern = 'SBS', barretts_final.sigs$Sigs)]
sigs.id <- barretts_final.sigs$Sigs[grepl(pattern = 'ID', barretts_final.sigs$Sigs)]

load("/Users/jao/Desktop/MSc_Project/Updated/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.sbs <- colnames(sigs.allSamples)


# Load hg19 signatures
signatures.sbs96_hg19 <- read.table('MutationalSignatureReferences/COSMIC_v3.3.1_SBS_GRCh37.txt', h=T)
sig_sbs96_types <- signatures.sbs96_hg19$Type
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, 2:ncol(signatures.sbs96_hg19)]
signatures.sbs96_hg19 <- data.frame(t(signatures.sbs96_hg19))
colnames(signatures.sbs96_hg19) <- sig_sbs96_types
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, colnames(signatures.cosmic)]

# Load COSMIC indel signatures
signatures.id83 <- read.table('MutationalSignatureReferences/COSMIC_v3.3_ID_GRCh37.txt', h=T)
sig_id83_types <- signatures.id83$Type
signatures.id83 <- signatures.id83[, 2:ncol(signatures.id83)]
signatures.id83 <- data.frame(t(signatures.id83))
colnames(signatures.id83) <- sig_id83_types

# Create id and sbs inputs
sigs.id83_input <- tcga_tally[,97:179]
sigs.id83_input <- sigs.id83_input[,colnames(signatures.id83)]
#sigs.id83_input <- sigs.id83_input[rowSums(sigs.id83_input) > 0, ]
sigs.id83_input <- sweep(sigs.id83_input, 1, rowSums(sigs.id83_input), FUN = "/")

sigs.sbs96_input <- tcga_tally[,1:96]
sigs.sbs96_input <- sigs.sbs96_input[,colnames(signatures.sbs96_hg19)]
#sigs.sbs96_input <- sigs.sbs96_input[rownames(sigs.sbs96_input) %in% rownames(sigs.id83_input), ] #85 samples

# ---- Run deconstructSigs

# Run deconstructSigs on each input tp calculate signature proportions
run_deconstructSigs <- function(sigs.input, sig_type = 'SBS') {
  
  if (sig_type == 'SBS') {
    sig_ref = signatures.sbs96_hg19[sigs.sbs, ] # only SBS signatures in ICGC
    print('Calculating SBS signature contributions...')
  } else if (sig_type == 'ID') {
    sig_ref = signatures.id83[sigs.id, ] # only ID signatures in ICGC
    print('Calculating ID signature contributions...')
  } else {
    print('Set sig_type to SBS or ID. Using SBS signatures...')
    sig_ref = signatures.sbs96_hg19[sigs.sbs, ]
  }
  
  sigs_out <- NULL
  for (i in 1:nrow(sigs.input)) {
    sample.i <- rownames(sigs.input)[i]
    print(paste0(sig_type, ' contributions, Sample ', i, ' of ', nrow(sigs.input), ': ', sample.i))
    
    sigs_i <- whichSignatures(
      tumor.ref = sigs.input,
      signatures.ref = sig_ref,
      sample.id = sample.i,
      contexts.needed = TRUE,
      signature.cutoff = 0,
      tri.counts.method = ifelse(sig_type == 'SBS', 'genome', 'default')
    )
    sigs_out <- rbind(sigs_out, sigs_i$weights)
  }
  
  return(sigs_out)
  
}
sigs.sbs96 <- run_deconstructSigs(sigs.sbs96_input, 'SBS')
sigs.id83 <- run_deconstructSigs(sigs.id83_input, 'ID')
tcga_sigs_complete <- cbind(sigs.sbs96, sigs.id83)


print("Saving")
save(tcga_sigs_complete,
     file = 'TCGA_deconstructSigs.Rdata')
