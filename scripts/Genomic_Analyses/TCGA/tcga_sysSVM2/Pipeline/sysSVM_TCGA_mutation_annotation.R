# Generate TCGA small somatic mutation annnotation sysSVM2 input - Server
source("/mnt/data/jao/sysSVM/sysSVM2/R/annotation_functions.R")

print("9-12-2024 run")

# ----- Small somatic mutation annotation - Server
vcf_files <- list.files(path = "/mnt/data/jao/sysSVM/sysSVM2/maf2vcf/vcf_tcga", pattern = "\\.vcf$", full.names = TRUE)
print(vcf_files)
ssms_annotated_list <- list()

for (i in seq_along(vcf_files)) {
  
  vcf <- vcf_files[i]
  
  sampleID <- sub("\\.vcf$", "", basename(vcf))
  
  ssms_annotated_list[[i]] <- annotate_ssms(
    vcf_fn = vcf, 
    sample = sampleID, 
    annovar_dir = "annovar", 
    genome_version = "hg38", 
    gene_aliases_entrez = "annotation_reference_files/gene_aliases_entrez.tsv", 
    hotspots = "annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv"
  )
  cat("VCF processed for sample:", sampleID, "\n")
}
ssms_annotated <- do.call(rbind, ssms_annotated_list)

save(ssms_annotated, file = "updated_tcga_ssms_annotated.RData")
