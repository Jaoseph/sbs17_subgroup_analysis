##### Processing of OCCAMS Barrett's samples only to infer SBS and indel mutation type contributions\

library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

#setwd("/Users/jao/Desktop/MSc_Project/Barretts_Classifier")

#----- Organise into MAF-readable format

# Single Nucleotide Variation
load("barretts.vcf.RData")
barretts.vcf$SampleID <- sub("_vs_.*", "", barretts.vcf$Sample) 
barretts.vcf$`#CHROM` <- gsub("chr", "", barretts.vcf$`#CHROM`)

# ONLY EXTRACT BARRETT'S SAMPLES

load("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/BarrettsNaiveComparison/annotation.sampleIDs.RData")
BE_samples <- annotation.sampleIDs$TumourID[annotation.sampleIDs$Category == "Barretts"]
barretts.vcf <- barretts.vcf[barretts.vcf$SampleID %in% BE_samples,]
length(unique(barretts.vcf$SampleID)) #check should be 127 samples


barretts_snv <- data.frame(
  Tumor_Sample_Barcode = barretts.vcf$SampleID,
  Chromosome = barretts.vcf$`#CHROM`,
  Start_Position = barretts.vcf$POS,
  End_Position = barretts.vcf$POS,
  Variant_Type = "SNP",
  Reference_Allele = barretts.vcf$REF,
  Tumor_Seq_Allele2 = barretts.vcf$ALT)

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
  Chromosome = oac_indels$Chromosome,
  Start_Position = oac_indels$Position,
  End_Position = oac_indels$Position,
  Variant_Type = sapply(oac_indels$REF,
                        function(x) ifelse(x == "-", "INS", "DEL")),
  Reference_Allele = oac_indels$REF,
  Tumor_Seq_Allele2 = oac_indels$ALT)

barretts_indels <- barretts_indels[barretts_indels$Tumor_Sample_Barcode %in% unique(barretts_snv$Tumor_Sample_Barcode),]

barretts_input <- rbind(barretts_snv, barretts_indels)
length(unique(barretts_input$Tumor_Sample_Barcode)) # 127 Samples

# Generate MAF
barretts_maf <- read_maf_minimal(barretts_input)

### Run mt_tally() from sigminer package to collate mutation type contributions

mt_tally.barretts <- sig_tally(
  object = barretts_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  useSyn = TRUE
)

### Saving 

print("Saving mt_tally.barretts")
save(mt_tally.barretts, 
     file = 'mt_tally.barretts.Rdata')


