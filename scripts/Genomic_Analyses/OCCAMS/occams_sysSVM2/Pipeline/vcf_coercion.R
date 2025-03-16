# ----- VCF_coercion (vcf file per sample) for sysSVM2 input - Server

# Load Sig17 annotations 
load("OCCAMS_CombinedresultsSummary.Rdata")

# ----- Small somatic mutation annotation
load("mutational_OCCAMS.Rdata")
mutational_OCCAMS <- mutational_OCCAMS[mutational_OCCAMS$SampleID %in% occams_results_df$Sample,] # Only extract samples with Sig17 annotations
mutational_OCCAMS <- mutational_OCCAMS[!mutational_OCCAMS$SampleID %in% c("SLX-18929_UDP0021", "SLX-18929_UDP0063"),] #Remove CNV missing samples
length(unique(mutational_OCCAMS$SampleID)) # Check # samples : 719 Samples

output_directory <- "vcf_files"

OCCAMS_samples <- unique(mutational_OCCAMS$SampleID)

for (sample in OCCAMS_samples) {
  
  # extract only the sample in current iteration
  sample_mutations <- mutational_OCCAMS[mutational_OCCAMS$SampleID == sample,]
  
  # Coerce into VCF format
  vcf_data <- data.frame(
    CHROM = sample_mutations$chr,       # Chromosome
    POS = sample_mutations$pos,         # Position
    ID = ".",                            # No variant ID available
    REF = sample_mutations$ref,         # Reference base
    ALT = sample_mutations$mut,         # Alternate base
    QUAL = ".",                          # Placeholder for quality score
    FILTER = ".",                        # Placeholder for filter
    INFO = "."                           # Placeholder for additional info
  )
  
  # Write VCF header
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    paste0("##source=mutational_OCCAMS"),
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )
  
  # Save as .vcf
  vcf_file <- file.path(output_directory, paste0(sample, ".vcf"))
  writeLines(vcf_header, vcf_file)  # Write the header
  write.table(
    vcf_data, file = vcf_file, sep = "\t", quote = FALSE,
    col.names = FALSE, row.names = FALSE, append = TRUE
  )
  cat("VCF file saved for sample:", sample, "\n")
}
