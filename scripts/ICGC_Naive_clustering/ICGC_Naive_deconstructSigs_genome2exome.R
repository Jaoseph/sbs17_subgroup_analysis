
##### Contribution of ICGC WGS-associated signatures in ICGC using deconstructSigs
setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")


# Load libraries
library(deconstructSigs)

# Load ICGC signature contributions
#   Only signatures appearing in >1% samples will be included
icgc_final.sigs <- read.table('NaiveOnly_sigProfilerCont.txt', header = TRUE)
icgc_final.sigs <- icgc_final.sigs[icgc_final.sigs$Prop_present > 0.01, ]
sigs.sbs <- icgc_final.sigs$Sigs[grepl(pattern = 'SBS', icgc_final.sigs$Sigs)]
sigs.id <- icgc_final.sigs$Sigs[grepl(pattern = 'ID', icgc_final.sigs$Sigs)]
load("/Users/jao/Desktop/MSc_Project/Updated/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.allSamples <- colnames(sigs.allSamples)
sigs.sbs <- unique(c(sigs.sbs, sigs.allSamples))

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
load('ICGC_Naive_mt_tally.Rdata')
sigs.sbs96_input <- as.data.frame(mt_tally.icgc_naive$SBS_96)
sigs.sbs96_input <- sigs.sbs96_input[,colnames(signatures.sbs96_hg19)]

sigs.id83_input <- as.data.frame(mt_tally.icgc_naive$ID_83)
sigs.id83_input <- sigs.id83_input[,colnames(signatures.id83)]

#remove samples with no INDELs
no_ID_samples <- read.table('no_ID_samples.txt', h=T)
no_ID_samples <- no_ID_samples$x

sigs.sbs96_input <- sigs.sbs96_input[!rownames(sigs.sbs96_input) %in% no_ID_samples,]
sigs.id83_input <- sigs.id83_input[!rownames(sigs.id83_input) %in% no_ID_samples,]

# Normalise ID83 counts
icgc_idCounts <- read.table('ICGC_Naive_indelCounts.txt')
icgc_idCounts$genome_to_exome <- icgc_idCounts$exome/icgc_idCounts$wgs
icgc_idCounts <- icgc_idCounts[colnames(sigs.id83_input), ]

sigs.id83_input <- as.data.frame(t(apply(sigs.id83_input, 1, function(x) x * icgc_idCounts$genome_to_exome)))

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
      tri.counts.method = ifelse(sig_type == 'SBS', 'genome2exome', 'default')
    )
    sigs_out <- rbind(sigs_out, sigs_i$weights)
  }
  
  return(sigs_out)
  
}

sigs.id83 <- run_deconstructSigs(sigs.id83_input, 'ID')
sigs.sbs96 <- run_deconstructSigs(sigs.sbs96_input, 'SBS')
sigs_complete <- cbind(sigs.sbs96, sigs.id83)

print("Saving")
save(sigs_complete,
     file = 'ICGC_Naive_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')
