##### Processing of Barretts samples only to infer SBS and indel mutation type contributions

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

print("Running Barretts_IndelCounts.R, 23_OCT")

# Set working directory
#setwd("/Users/jao/Desktop/MSc_Project/Barretts_Classifier")

# Indels

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

barretts_indels <- data.frame(
  Tumor_Sample_Barcode = oac_indels$SampleID,
  Hugo_Symbol = NA,
  Chromosome = oac_indels$Chromosome,
  Start_Position = oac_indels$Position,
  End_Position = oac_indels$Position,
  Variant_Classification = oac_indels$V7,
  Variant_Type = sapply(oac_indels$REF,
                        function(x) ifelse(x == "-", "INS", "DEL")),
  Reference_Allele = oac_indels$REF,
  Tumor_Seq_Allele2 = oac_indels$ALT)


# ONLY EXTRACT BARRETT'S SAMPLES
load("annotation.sampleIDs.RData")
BE_samples <- annotation.sampleIDs$TumourID[annotation.sampleIDs$Category == "Barretts"]
barretts_indels <- barretts_indels[barretts_indels$Tumor_Sample_Barcode %in% BE_samples,]
length(unique(barretts_indels$Tumor_Sample_Barcode)) #check should be 127 samples

# Generate maf for whole-genome

barretts_maf <- read.maf(maf = barretts_indels,
                         vc_nonSyn = names(table(barretts_indels$Variant_Classification)))

syn_indel <- c("inframe_insertion","inframe_deletion",
               "frameshift_variant", "frameshift_variant,NMD_transcript_variant",
               "TF_binding_site_variant,TFBS_ablation",
               "frameshift_variant,splice_region_variant",
               "stop_gained,frameshift_variant",
               "stop_gained,frameshift_variant,NMD_transcript_variant",
               "frameshift_variant,splice_region_variant,NMD_transcript_variant",
               "inframe_insertion,NMD_transcript_variant",
               "splice_acceptor_variant,frameshift_variant",
               "frameshift_variant,initiator_codon_variant",
               "inframe_deletion,incomplete_terminal_codon_variant",
               "inframe_deletion,NMD_transcript_variant",
               "frameshift_variant,stop_lost,NMD_transcript_variant",
               "inframe_deletion,splice_region_variant",
               "frameshift_variant,stop_lost",
               "incomplete_terminal_codon_variant,coding_sequence_variant,3_prime_UTR_variant",
               "inframe_deletion,splice_region_variant,NMD_transcript_variant",
               "frameshift_variant,stop_lost,splice_region_variant,NMD_transcript_variant",
               "stop_gained,frameshift_variant,splice_region_variant",
               "inframe_insertion,splice_region_variant",
               "frameshift_variant,stop_retained_variant",
               "stop_gained,inframe_insertion",
               "frameshift_variant,stop_lost,splice_region_variant",
               "inframe_insertion,incomplete_terminal_codon_variant,coding_sequence_variant",
               "frameshift_variant,stop_retained_variant,NMD_transcript_variant",
               "inframe_insertion,splice_region_variant,NMD_transcript_variant",
               "stop_gained,inframe_deletion",
               "initiator_codon_variant,inframe_deletion",
               "stop_gained,inframe_insertion,NMD_transcript_variant",
               "stop_gained,inframe_deletion,NMD_transcript_variant")

# Generate maf for exome

barretts_exome_indels <- barretts_indels[barretts_indels$Variant_Classification %in% syn_indel,]

barretts_exome_maf <- read.maf(maf = barretts_exome_indels,
                               vc_nonSyn = names(table(barretts_exome_indels$Variant_Classification)))


# Calculate ID-83 counts for each MAF file

print("Running ID mt_tally() for wgs")
mt_tally.barretts_wgs_id <- sig_tally(
  object = barretts_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

print("Running ID mt_tally() for downsampled wxs")
mt_tally.barretts_exome_id <- sig_tally(
  object = barretts_exome_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

# Count total mutations of each ID-83 type
barretts_indelCounts <- data.frame(
  wgs = apply(mt_tally.barretts_wgs_id$all_matrices$ID_83, 2, sum),
  exome = apply(mt_tally.barretts_exome_id$all_matrices$ID_83, 2, sum)
)

print("writing Barretts_indelCounts.txt")
write.table(barretts_indelCounts, file = 'Barretts_indelCounts.txt')
