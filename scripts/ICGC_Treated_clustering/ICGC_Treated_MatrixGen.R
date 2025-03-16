##### Processing of ICGC Treated samples only to infer SBS and indel mutation type contributions

library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

# Set working directory
#setwd("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only")

#----- Extract Chemotherapy Naive Samples data

# Read chemo-treated samples
#treated_annotations <- read.csv("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/ICGC_Treated_Samples.csv")
treated_annotations <- read.csv("ICGC_Treated_Samples.csv")

# Identify projects
#projects <- unique(treated_annotations$project_code) ; projects
#file_path <- "/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/ICGC_Treated_Samples/"
datasets <- list.files(pattern = "simple_somatic_mutation.open.*.tsv.gz")

# Extract Somatic Mutations
print("Starting Loop")

icgc_treated <- data.frame()

for (i in 1:length(datasets)){
  file <- datasets[i]
  data <- read.table(file, sep = '\t', header = TRUE)
  data_wgs <- data[data$sequencing_strategy == 'WGS', ]
  data_wgs <- data_wgs[!duplicated(data_wgs$icgc_mutation_id), ]
  icgc_treated <- rbind(icgc_treated, data_wgs)
  
  #Output reporting
  project_identifier <- gsub("simple_somatic_mutation.open.|.tsv.gz", "", file)
  print(paste0(i,' - Extracted WGS-specific somatic mutations for ', project_identifier))
}

# Extract ChemoNaive samples only
treated_samples <- icgc_treated[icgc_treated$icgc_donor_id %in% treated_annotations$icgc_donor_id,]

#----- Organise into MAF-readable format

print("Generating icgc_treated_input")

icgc_treated_input <- data.frame(
  Tumor_Sample_Barcode = treated_samples$icgc_donor_id,
  Hugo_Symbol = NA,
  Chromosome = treated_samples$chromosome,
  Start_position = treated_samples$chromosome_start,
  End_position = treated_samples$chromosome_end,
  Variant_Classification = treated_samples$consequence_type,
  Variant_Type = sapply(treated_samples$mutation_type,
                        function(x) ifelse(x == 'single base substitution', 'SNP',
                                           ifelse(x == 'insertion of <=200bp', 'INS',
                                                  ifelse(x == 'deletion of <=200bp', 'DEL',NA)))),
  Reference_Allele = treated_samples$reference_genome_allele,
  Tumor_Seq_Allele2 = treated_samples$mutated_to_allele
)

icgc_treated_maf <- read.maf(maf = icgc_treated_input,
                           vc_nonSyn = names(table(icgc_treated_input$Variant_Classification)))


# Run mt_tally() from sigminer package to collate mutation type contributions

print("Running mt_tally")

mt_tally.icgc_treated <- sig_tally(
  object = icgc_treated_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  useSyn = TRUE
)

save(mt_tally.icgc_treated, 
     file = 'icgc_treated_mt_tally.Rdata')