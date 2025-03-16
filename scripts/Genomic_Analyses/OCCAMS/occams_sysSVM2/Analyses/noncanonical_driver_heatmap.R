setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
library(dplyr)
library(wesanderson)
library(ComplexHeatmap)
library(tidyverse)
library(openxlsx)

# Load driver predictions
load("sysSVM2_input.Rdata")
sysSVM2_input <- subset(sysSVM2_input, select = c("sample", "entrez", "no_ALL_muts", "no_TRUNC_muts",
                                                  "no_NTDam_muts", "no_GOF_muts", "Copy_number", "CNVGain", "CNVLoss"))

# A list of predicted drivers in each sample
drivers_list <- readRDS("test_sysSVM2/drivers_toppedUp.rds")

drivers_list <- left_join(drivers_list, sysSVM2_input, by = c("sample", "entrez"))

# Load Sig17 annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample, Phenotype_Assigned))
names(occams_results_df) <- c("sample", "Phenotype_Assigned")

drivers_list <- left_join(drivers_list, occams_results_df, by = "sample")


heatmap <- drivers_list %>%
  group_by(Phenotype_Assigned) %>%
  mutate(n_samples = n_distinct(sample))%>%
  filter(canonical_driver == FALSE) %>%
  group_by(Phenotype_Assigned, symbol) %>%
  summarise(n_samples = n())

occams_canonicaldrivers <- data.frame(
  Phenotype_Assigned = heatmap$Phenotype_Assigned,
  Driver = heatmap$symbol)
write.xlsx(occams_canonicaldrivers, "Analysis/OCCAMS_noncanonicaldrivers_all.xlsx", rowNames = FALSE)

heatmap <- heatmap[heatmap$n_samples >= 2,] #for non-canonical, only include drivers found in at least 2 samples

occams_canonicaldrivers <- data.frame(
  Phenotype_Assigned = heatmap$Phenotype_Assigned,
  Driver = heatmap$symbol)
write.xlsx(occams_canonicaldrivers, "Analysis/OCCAMS_noncanonicaldrivers_2samples.xlsx", rowNames = FALSE)

heatmap2 <- heatmap %>%
  pivot_wider(names_from = Phenotype_Assigned, values_from = n_samples) %>%
  column_to_rownames("symbol")
heatmap2[is.na(heatmap2)] <- 0
heatmap2 <- heatmap2[,c(1,2,4,3)]

# Generate Heatmap Plot
group_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

phenotype_groups <- factor(colnames(heatmap2), 
                           levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-"))

#column annotations
col_ann <- HeatmapAnnotation(
  group = phenotype_groups,  
  col = list(group = group_colors),
  annotation_legend_param = list(
    title = "Phenotype"
  ), show_annotation_name = FALSE)

#row annotations
can_driver <- data.frame(
  gene = rownames(heatmap2)
)
esophageal <- read.delim("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/OACdrivers_Frankell.txt", header = FALSE)$V1 
can_driver$eac <- ifelse(can_driver$gene %in% esophageal, "yes", "no")

binary_colors <- c("yes" = "black", "no" = "white")

#row_ann <- rowAnnotation(
#  group = can_driver$eac, 
#  col = list(group = binary_colors))

row_ann <- rowAnnotation(
  eso = anno_simple(can_driver$eac, col = binary_colors), show_annotation_name = FALSE)

# Define the heatmap
ht <- Heatmap(
  as.matrix(heatmap2),
  col = circlize::colorRamp2(c(0, 5), c("white", "red")),
  heatmap_legend_param = list(title = "Proportion"),
  show_heatmap_legend = FALSE,
  cluster_rows = TRUE,
  column_split = phenotype_groups,
  column_title_gp = grid::gpar(fontsize = 0),
  row_names_gp = grid::gpar(fontsize = 5, fontface = "bold"),
  show_column_names = FALSE,
  top_annotation = col_ann,
  right_annotation = row_ann,
  cluster_columns = FALSE
)

# Define custom legends
ht_legend <- Legend(title = "Number of Samples", col_fun = circlize::colorRamp2(c(0, 10), c("white", "red")))
anno_legends <- list(
  Legend(title = "Previous EAC Driver", at = c("yes", "no"), labels = c("Yes", "No"), legend_gp = gpar(fill = c("black", "white")))
  #,Legend(title = "Phenotype", at = phenotype_groups, legend_gp = gpar(fill = group_colors))
)

# Combine all legends into a single object
combined_legends <- packLegend(
  ht_legend,
  anno_legends[[1]],
  #anno_legends[[2]],
  direction = "vertical"
)

pdf("Analysis/NonCanonical_Driver_Heatmap.pdf", h = 10, w = 7)
draw(
  ht,
  annotation_legend_list = combined_legends,
  annotation_legend_side = "right",
  heatmap_legend_side = "right"
)
dev.off()


# ----- Generate Spreadsheet

supp_table <- heatmap2
supp_table$gene <- rownames(supp_table)

supp_table <- merge(supp_table, can_driver, by = "gene")
supp_table <- supp_table %>%
  mutate(priority = ifelse(`NaiveLike.Sig17+` > 0 & `BarrettsLike.Sig17+` == 0 & `TreatedLike.Sig17+` == 0 & `Sig17-` == 0, 1, 2)) %>%
  arrange(priority, desc(`NaiveLike.Sig17+`), desc(`Sig17-`)) %>% 
  select(-priority) 

# Save
library(openxlsx)
write.xlsx(supp_table, "occams_noncanHeatmap_supp_table.xlsx", sheetName = "Table S", rowNames=TRUE)

