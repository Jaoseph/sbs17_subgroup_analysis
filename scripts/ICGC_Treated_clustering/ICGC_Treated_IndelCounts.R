##### Processing of ICGC Treated samples only to infer SBS and indel mutation type contributions

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)


# Set working directory
#setwd("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only")

#----- Extract Chemotherapy Treated Samples data

# Read Chemo Treated samples
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

# Extract Chemo Treated samples only
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

# Remove mutations in intronic regions

icgc_exome_input <- icgc_treated_input[icgc_treated_input$Variant_Classification != 'intergenic_region',]

icgc_maf <- read.maf(maf = icgc_treated_input,
                     vc_nonSyn = names(table(icgc_treated_input$Variant_Classification)))

icgc_exome_maf <- read.maf(maf = icgc_exome_input,
                           vc_nonSyn = names(table(icgc_exome_input$Variant_Classification)))


# Calculate ID-83 counts for each MAF file

print("Running ID mt_tally() for wgs")
mt_tally.icgc_wgs_id <- sig_tally(
  object = icgc_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

print("Running ID mt_tally() for downsampled wxs")
mt_tally.icgc_exome_id <- sig_tally(
  object = icgc_exome_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

# Count total mutations of each ID-83 type
icgc_indelCounts <- data.frame(
  wgs = apply(mt_tally.icgc_wgs_id$all_matrices$ID_83, 2, sum),
  exome = apply(mt_tally.icgc_exome_id$all_matrices$ID_83, 2, sum)
)


print("writing ICGC_indelCounts.txt")
write.table(icgc_indelCounts, file = 'icgc_treated_indelCounts.txt')