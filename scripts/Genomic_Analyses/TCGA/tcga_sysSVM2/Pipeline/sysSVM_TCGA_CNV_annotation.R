setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/")
source("R/annotation_functions.R")

# ---- CNV Segment annotations; TCGA

# Read ASCAT TCGA-ESCA CNV Seg data
seg_data <- read.table("TCGA/esca_tcga_segments.seg", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cnv_segments <- seg_data %>%
  mutate(
    chromosome = chrom,
    start = loc.start,
    end = loc.end,
    copy_number = round(2^(seg.mean) * 2) # Convert seg.mean to copy number
  ) %>%
  select(sample = ID, chromosome, start, end, copy_number)

# Extract only EAC samples
tcga_annotation <- read.delim('TCGA/esca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep="\t")
tcga_annotation <- tcga_annotation[tcga_annotation$Cancer.Type.Detailed == "Esophageal Adenocarcinoma",]
tcga_eac_ID_index <- tcga_annotation$Patient.ID

cnv_segments$sample <- sapply(cnv_segments$sample,
                                 function(x) substr(x,1,12))
cnv_segments <- cnv_segments[cnv_segments$sample %in% tcga_eac_ID_index,]

# generate sysSVM2 CNV inout
cnvs_annotated = annotate_cnvs(
  cnv_segments, 
  bedtools_bin_dir = "/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/bedtools2/bin",
  gene_coords = "annotation_reference_files/gene_coords_hg19.tsv"
)

