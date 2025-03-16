setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/Analysis")
library(dplyr)
library(wesanderson)

# Load driver predictions
predictions <- readRDS("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/drivers_toppedUp.rds")
predictions <- subset(predictions, select = c(sample, canonical_driver, symbol))

# Load Sig17 annotations for tcga smaples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga$sample <- rownames(ann_tcga)
ann_tcga <- subset(ann_tcga, select = c(sample, Phenotype_Assigned))
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",] #remove BElike sample.

driver_df <- left_join(predictions, ann_tcga, by = "sample")
driver_df <- driver_df[, c(4,1,3,2)] ; driver_df <- driver_df[driver_df$sample != "TCGA-L5-A4OT",] #remove BElike sample.

library(ggvenn)
library(tidyverse)
library(wesanderson)

venn_df <- split(driver_df$symbol, driver_df$Phenotype_Assigned)
save(venn_df, file = "tcga_drivers_list.Rdata")

pdf("TCGA_Drivers_Venn.pdf")
ggvenn(venn_df, 
       fill_color = c(wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[2]),
       stroke_size = 0.5, set_name_size = 4, text_size = 3) 
dev.off()



# ----- Identifying which genes are potential drivers in which subgroups

# Load driver predictions
predictions <- readRDS("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2/TCGA/drivers_toppedUp.rds")
predictions <- subset(predictions, select = c(sample, canonical_driver, symbol))

# Load Sig17 annotations for tcga smaples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga$sample <- rownames(ann_tcga)
ann_tcga <- subset(ann_tcga, select = c(sample, Phenotype_Assigned))
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",] #remove BElike sample. (TCGA-L5-A4OT)

driver_df <- left_join(predictions, ann_tcga, by = "sample")
driver_df <- driver_df[, c(4,1,3,2)] ; driver_df <- driver_df[driver_df$sample != "TCGA-L5-A4OT",]

# Generate Matrix
matrix <- array(FALSE,c(length(unique(driver_df$symbol)),length(unique(driver_df$Phenotype_Assigned))))
colnames(matrix) <- unique(driver_df$Phenotype_Assigned)
rownames(matrix) <- unique(driver_df$symbol)

for(i in 1:nrow(driver_df)){
  matrix[driver_df[i,]$symbol,driver_df[i,]$Phenotype_Assigned] <- TRUE
}
matrix <- as.data.frame(matrix)

matrix$n_groups <- rowSums(matrix)
matrix$category <- apply(matrix[,c(1:2)], 1, function(x) paste(colnames(matrix[,c(1:2)])[x == TRUE], collapse = ":"))
matrix <- matrix %>%
  arrange(category)

# Save
library(openxlsx)
write.xlsx(matrix, "Gene_Presence_TCGA.xlsx", sheetName = "Table S", rowNames=TRUE)
