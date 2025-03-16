setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/Analysis")
library(dplyr)
library(wesanderson)
library(ComplexHeatmap)
library(tidyverse)
library(openxlsx)

# Load annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/sysSVM2_tcga_input.Rdata")

# A list of predicted drivers in each sample
drivers_list <- readRDS("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/drivers_toppedUp.rds")
drivers_list <- left_join(drivers_list, molecular_data, by = c("sample", "entrez"))

# Load Sig17 annotations for tcga smaples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga$sample <- rownames(ann_tcga)
ann_tcga <- subset(ann_tcga, select = c(sample, Phenotype_Assigned))
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",] #remove BElike sample.

drivers_list <- left_join(drivers_list, ann_tcga, by = "sample")
drivers_list <- drivers_list[drivers_list$sample != "TCGA-L5-A4OT",] #remove BElike sample.

# Format for Heatmap plotting

heatmap <- drivers_list %>%
  group_by(Phenotype_Assigned) %>%
  mutate(n_samples = n_distinct(sample))%>%
  filter(canonical_driver == FALSE) %>%
  group_by(Phenotype_Assigned, symbol) %>%
  summarise(n_samples = n())

tcga_noncanonicaldrivers <- data.frame(
  Phenotype_Assigned = heatmap$Phenotype_Assigned,
  Driver = heatmap$symbol)
write.xlsx(tcga_noncanonicaldrivers, "tcga_noncanonicaldrivers_allsamples.xlsx", rowNames = FALSE)

#heatmap <- heatmap[heatmap$n_samples >= 2,] #for non-canonical, majority is single.

heatmap2 <- heatmap %>%
  pivot_wider(names_from = Phenotype_Assigned, values_from = n_samples) %>%
  column_to_rownames("symbol")
heatmap2[is.na(heatmap2)] <- 0
write.xlsx(heatmap2, "tcga_noncanonicaldrivers.xlsx", rowNames = TRUE)

# ---- Generate Heatmap Plot

group_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4])
phenotype_groups <- factor(colnames(heatmap2), 
                           levels = c("NaiveLike.Sig17+", "Sig17-"))

#column annotations
col_ann <- HeatmapAnnotation(
  group = phenotype_groups,  
  col = list(group = group_colors),
  annotation_legend_param = list(
    title = "Phenotype"
  ), show_annotation_name = FALSE)

#row annotations
can_driver <- data.frame(
  gene = rownames(heatmap2))
esophageal <- read.delim("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/OACdrivers_Frankell.txt", header = FALSE)$V1 
can_driver$eac <- ifelse(can_driver$gene %in% esophageal, "yes", "no")

binary_colors <- c("yes" = "black", "no" = "white")

row_ann <- rowAnnotation(
  eso = anno_simple(can_driver$eac, col = binary_colors), show_annotation_name = FALSE)

# Define the heatmap
ht <- Heatmap(
  as.matrix(heatmap2),
  col = circlize::colorRamp2(c(0, 5), c("white", "red")),
  heatmap_legend_param = list(title = "Proportion"),
  show_heatmap_legend = FALSE,
  cluster_rows = TRUE,
  column_split = 4,
  row_names_gp = grid::gpar(fontsize = 6, fontface = "bold"),
  show_column_names = FALSE,
  top_annotation = col_ann,
  right_annotation = row_ann)

# Define custom legends
ht_legend <- Legend(title = "Number of Samples", col_fun = circlize::colorRamp2(c(0, 3), c("white", "red")))
anno_legends <- list(
  Legend(title = "Previous EAC Driver", at = c("yes", "no"), labels = c("Yes", "No"), legend_gp = gpar(fill = c("black", "white")))
  #,Legend(title = "Phenotype", at = phenotype_groups, legend_gp = gpar(fill = group_colors))
)

# Combine all legends into a single object
combined_legends <- packLegend(
  ht_legend,
  anno_legends[[1]],
  #anno_legends[[2]],
  direction = "vertical")


#pdf("TCGA_NonCanonical_Driver_Heatmap.pdf", h = 10, w = 7)
draw(ht, annotation_legend_list = combined_legends,
     annotation_legend_side = "right",
     heatmap_legend_side = "right")
#dev.off()


# ----- Generate Spreadsheet

supp_table <- heatmap2
supp_table$gene <- rownames(supp_table)

supp_table <- merge(supp_table, can_driver, by = "gene")
supp_table <- supp_table %>%
  arrange(desc(`NaiveLike.Sig17+`), `Sig17-`)

# Save
library(openxlsx)
write.xlsx(supp_table, "tcga_noncanHeatmap_supp_table.xlsx", sheetName = "Table S", rowNames=TRUE)
