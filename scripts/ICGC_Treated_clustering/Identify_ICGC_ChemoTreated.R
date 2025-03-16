##### Identify Chemo Treated Samples

library(dplyr)
library(purrr)
library(readr)
library(stringr)

#set working directory
setwd("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only")

#Identify all files
file_path <- "/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/Specimen/all_projects"
tsv_files <- list.files(file_path, pattern = "\\.tsv$", full.names = TRUE)

#Extract all samples with treatment information
combined_data <- tsv_files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE, col_types = cols(.default = "c")), .id = "file_name")

combined_data <- combined_data[,-1]

# Add colnames  
header <- read_tsv("/Users/jao/Desktop/MSc_Project/ChemoTreated_Only/Specimen/specimen.tsv.gz")
colnames(combined_data) <- colnames(header)

combined_data <- combined_data %>%
  filter(!str_detect(specimen_type, "Normal") & !str_detect(specimen_type, "Metastatic"))

# Extract Chemo Naive 

chemo_treatment <- c("chemotherapy", "combined chemo+radiation therapy", "combined chemo+immunotherapy")

treated_samples <- combined_data[combined_data$specimen_donor_treatment_type %in% chemo_treatment,]
treated_samples <- treated_samples[!duplicated(treated_samples$icgc_donor_id),]
treated_samples <- subset(treated_samples, !is.na(specimen_donor_treatment_type)) #374 samples in total across

write.csv(treated_samples, file = "ICGC_Treated_Samples.csv", row.names = FALSE)

