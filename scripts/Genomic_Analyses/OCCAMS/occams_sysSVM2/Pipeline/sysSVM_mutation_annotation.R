# Generate small somatic mutation annnotation sysSVM2 input - Server
source("R/annotation_functions.R")

# ----- Small somatic mutation annotation - Server
vcf_files <- list.files(path = "vcf_files", pattern = "\\.vcf$", full.names = TRUE)

ssms_annotated_list <- list()

for (i in seq_along(vcf_files)) {
  
  vcf <- vcf_files[i]
  
  sampleID <- sub("\\.vcf$", "", basename(vcf))
  
  ssms_annotated_list[[i]] <- annotate_ssms(
    vcf_fn = vcf, 
    sample = sampleID, 
    annovar_dir = "annovar", 
    genome_version = "hg19", 
    gene_aliases_entrez = "annotation_reference_files/gene_aliases_entrez.tsv", 
    hotspots = "annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv"
  )
  cat("VCF processed for sample:", sampleID, "\n")
}
ssms_annotated <- do.call(rbind, ssms_annotated_list)

save(ssms_annotated, file = "ssms_annotated.RData")
