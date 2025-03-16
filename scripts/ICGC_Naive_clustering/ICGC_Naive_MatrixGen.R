##### Processing of ICGC Naive samples only to infer SBS and indel mutation type contributions

library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Set working directory
#setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_ChemoNaive_Samples")

#----- Extract Chemotherapy Naive Samples data

# Read chemonaive samples
#naive_annotations <- read.csv("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_NaiveOnly_Samples.csv")
naive_annotations <- read.csv("ICGC_NaiveOnly_Samples.csv")

# Identify projects
#projects <- unique(naive_annotations$project_code) ; projects
#file_path <- "/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_ChemoNaive_Samples/"
datasets <- list.files(pattern = "simple_somatic_mutation.open.*.tsv.gz")

# Extract Somatic Mutations
print("Starting Loop")

icgc_naive <- data.frame()

for (i in 1:length(datasets)){
  file <- datasets[i]
  data <- read.table(file, sep = '\t', header = TRUE)
  data_wgs <- data[data$sequencing_strategy == 'WGS', ]
  data_wgs <- data_wgs[!duplicated(data_wgs$icgc_mutation_id), ]
  icgc_naive <- rbind(icgc_naive, data_wgs)
  
  #Output reporting
  project_identifier <- gsub("simple_somatic_mutation.open.|.tsv.gz", "", file)
  print(paste0(i,' - Extracted WGS-specific somatic mutations for ', project_identifier))
}

# Extract ChemoNaive samples only
naive_samples <- icgc_naive[icgc_naive$icgc_donor_id %in% naive_annotations$icgc_donor_id,]

#----- Organise into MAF-readable format

print("Generating icgc_naive_input")

icgc_naive_input <- data.frame(
  Tumor_Sample_Barcode = naive_samples$icgc_donor_id,
  Hugo_Symbol = NA,
  Chromosome = naive_samples$chromosome,
  Start_position = naive_samples$chromosome_start,
  End_position = naive_samples$chromosome_end,
  Variant_Classification = naive_samples$consequence_type,
  Variant_Type = sapply(naive_samples$mutation_type,
                        function(x) ifelse(x == 'single base substitution', 'SNP',
                                           ifelse(x == 'insertion of <=200bp', 'INS',
                                                  ifelse(x == 'deletion of <=200bp', 'DEL',NA)))),
  Reference_Allele = naive_samples$reference_genome_allele,
  Tumor_Seq_Allele2 = naive_samples$mutated_to_allele
)

icgc_naive_maf <- read.maf(maf = icgc_naive_input,
                     vc_nonSyn = names(table(icgc_naive_input$Variant_Classification)))


# Run mt_tally() from sigminer package to collate mutation type contributions

print("Running mt_tally")

mt_tally.icgc_naive <- sig_tally(
  object = icgc_naive_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  useSyn = TRUE
)

save(mt_tally.icgc_naive, 
     file = 'ICGC_Naive_mt_tally.Rdata')