setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
library(dplyr)

# Load driver predictions
predictions <- readRDS("test_sysSVM2/drivers_toppedUp.rds")
predictions <- subset(predictions, select = c(sample, canonical_driver, symbol))

# Load Sig17 annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample, Phenotype_Assigned))
names(occams_results_df) <- c("sample", "Phenotype_Assigned")

driver_df <- left_join(predictions, occams_results_df, by = "sample")
driver_df <- driver_df[, c(4,1,3,2)]
#names(driver_df) <- c("Signature Phenotype", "Sample", "Gene", "Canonical Driver")

library(ggvenn)
library(tidyverse)
library(wesanderson)

venn_df <- split(driver_df$symbol, driver_df$Phenotype_Assigned)
venn_df <- venn_df[c("BarrettsLike.Sig17+", "NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-")]
save(venn_df, file = "drivers_list.Rdata")

pdf("Analysis/Drivers_Venn.pdf",w = 7, h = 7)
ggvenn(venn_df, 
       fill_color = c(wes_palette("AsteroidCity3")[3], wes_palette("AsteroidCity3")[4], wes_palette("AsteroidCity3")[1], wes_palette("AsteroidCity3")[2]),
       stroke_size = 0.5, set_name_size = 4, text_size = 3.3) 
dev.off()



# ----- Identifying which genes are potential drivers in which subgroups

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
library(dplyr)

# Load driver predictions
predictions <- readRDS("test_sysSVM2/drivers_toppedUp.rds")
predictions <- subset(predictions, select = c(sample, canonical_driver, symbol))

# Load Sig17 annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample, Phenotype_Assigned))
names(occams_results_df) <- c("sample", "Phenotype_Assigned")

driver_df <- left_join(predictions, occams_results_df, by = "sample")
driver_df <- driver_df[, c(4,1,3,2)]

# Generate Matrix
matrix <- array(FALSE,c(length(unique(driver_df$symbol)),length(unique(driver_df$Phenotype_Assigned))))
colnames(matrix) <- unique(driver_df$Phenotype_Assigned)
rownames(matrix) <- unique(driver_df$symbol)

for(i in 1:nrow(driver_df)){
  matrix[driver_df[i,]$symbol,driver_df[i,]$Phenotype_Assigned] <- TRUE
}
matrix <- as.data.frame(matrix)

matrix$n_groups <- rowSums(matrix)
matrix$category <- apply(matrix[,c(1:4)], 1, function(x) paste(colnames(matrix[,c(1:4)])[x == TRUE], collapse = ":"))
matrix <- matrix %>%
  arrange(category)

# Save
library(openxlsx)
write.xlsx(matrix, "Analysis/Gene_Presence_4Groups.xlsx", sheetName = "Table S", rowNames=TRUE)
