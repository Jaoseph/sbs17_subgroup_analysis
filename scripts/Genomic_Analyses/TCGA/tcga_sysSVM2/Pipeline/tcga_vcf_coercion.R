# ----- VCF_coercion (vcf file per sample) for sysSVM2 input - Server

# Install TCGAbiolinks (Bioconductor)
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-ESCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
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

# Remove chr for Chromosome
tcga_mutations$Chromosome <- gsub("^chr", "", tcga_mutations$Chromosome)

# Loop to coerce.
output_directory <- "vcf_files_TCGA"

tcga_samples <- unique(tcga_mutations$Patient_ID)

for (sample in tcga_samples){
  # extract only the sample in current iteration
  sample_mutations <- tcga_mutations[tcga_mutations$Patient_ID == sample,]
  
  # Coerce into VCF format
  vcf_data <- data.frame(
    CHROM = sample_mutations$Chromosome,      
    POS = sample_mutations$Start_Position,        
    ID = ".",                           
    REF = sample_mutations$Reference_Allele,      
    ALT = sample_mutations$Tumor_Seq_Allele2,     
    QUAL = ".",                         
    FILTER = ".",                       
    INFO = ".")
  
  # write VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    paste0("##source=tcga_mutations"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
  
  # Save as .vcf
  vcf_file <- file.path(output_directory, paste0(sample, ".vcf"))
  writeLines(vcf_header, vcf_file)  # Write the header
  write.table(
    vcf_data, file = vcf_file, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE, append = TRUE
  )
  cat("VCF file saved for sample:", sample, "\n")
}