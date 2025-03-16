setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
source("R/annotation_functions.R")

load("OCCAMS_CombinedresultsSummary.Rdata")

# ----- CNV segment annotations

oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$sample <- sub("_vs_.*", "", oac_cna$sample)
oac_cna <- oac_cna[oac_cna$sample %in%  occams_results_df$Sample,]
# 719 samples only, missing 719 ("SLX-18929_UDP0021" "SLX-18929_UDP0063")

# Table with the following columns: sample; ploidy. Leave null if unavailable (assumes diploidy)
ploidy <- unique(oac_cna[,c("sample","ploidy")]) 

# Table with the following columns: sample; chromosome; start; end; and copy_number or segment_mean
oac_cna <- oac_cna[,c("sample", "chromosome", "start", "end", "total_cn")]
names(oac_cna) <- c("sample", "chromosome", "start", "end", "copy_number")

oac_cna$start <- as.numeric(oac_cna$start) 
oac_cna$end <- as.numeric(oac_cna$end) 
oac_cna$copy_number <- as.numeric(oac_cna$copy_number) 

cnvs_annotated = annotate_cnvs(
  oac_cna, 
  bedtools_bin_dir = "/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/bedtools2/bin",
  ploidy, 
  gene_coords = "annotation_reference_files/gene_coords_hg19.tsv"
)

save(cnvs_annotated, file = "cnvs_annotated.Rdata")
