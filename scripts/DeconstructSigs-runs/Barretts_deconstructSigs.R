setwd("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs")

# Load libraries
library(deconstructSigs)

# Load Barretts signature contributions
#   Only signatures appearing in >1% samples will be included
barretts_final.sigs <- read.table('EAC_sigProfilerCont.txt', header = TRUE)
barretts_final.sigs <- barretts_final.sigs[barretts_final.sigs$Prop_present > 0.01, ]
#sigs.sbs <- barretts_final.sigs$Sigs[grepl(pattern = 'SBS', barretts_final.sigs$Sigs)]
sigs.id <- barretts_final.sigs$Sigs[grepl(pattern = 'ID', barretts_final.sigs$Sigs)]

load("/Users/jao/Desktop/MSc_Project/Updated/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.sbs <- colnames(sigs.allSamples)

# Load hg19 signatures (for ICGC samples)
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

# Load in datasets and tweak into deconstructSigs inputs:
#   rows = samples, cols = signature contexts
#   order the same as signature data frames
load('/Users/jao/Desktop/MSc_Project/Barretts_Classifier/mt_tally.barretts.Rdata')
sigs.sbs96_input <- as.data.frame(mt_tally.barretts$SBS_96)
sigs.sbs96_input <- sigs.sbs96_input[,colnames(signatures.sbs96_hg19)]

sigs.id83_input <- as.data.frame(mt_tally.barretts$ID_83)
sigs.id83_input <- sigs.id83_input[,colnames(signatures.id83)]
sigs.id83_input <- sweep(sigs.id83_input, 1, rowSums(sigs.id83_input), FUN = "/")

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

sigs.id83 <- run_deconstructSigs(sigs.id83_input, 'ID')
sigs.sbs96 <- run_deconstructSigs(sigs.sbs96_input, 'SBS')
barretts_sigs_complete <- cbind(sigs.sbs96, sigs.id83)

print("Saving")
save(barretts_sigs_complete,
     file = 'Barretts_deconstructSigs.Rdata')
