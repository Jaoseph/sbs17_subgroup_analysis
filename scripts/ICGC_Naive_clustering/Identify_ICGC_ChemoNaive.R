##### Identify Chemo Naive Samples

library(dplyr)
library(purrr)
library(readr)
library(stringr)

#set working directory
setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")

#Identify all files
file_path <- "/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/Specimen/all_projects"
tsv_files <- list.files(file_path, pattern = "\\.tsv$", full.names = TRUE)

#Extract all samples with treatment information
combined_data <- tsv_files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE, col_types = cols(.default = "c")), .id = "file_name")

combined_data <- combined_data[,-1]

# Add colnames  
header <- read_tsv("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/Specimen/specimen.tsv.gz")
colnames(combined_data) <- colnames(header)

combined_data <- combined_data %>%
  filter(!str_detect(specimen_type, "Normal") & !str_detect(specimen_type, "Metastatic"))

# Extract Chemo Naive 

chemo_treatment <- c("chemotherapy", "combined chemo+radiation therapy", "combined chemo+immunotherapy")

naive_samples <- combined_data[!combined_data$specimen_donor_treatment_type %in% chemo_treatment,]
naive_samples <- naive_samples[!duplicated(naive_samples$icgc_donor_id),]
naive_samples <- subset(naive_samples, !is.na(specimen_donor_treatment_type)) #4844 samples in total across

write.csv(naive_samples, file = "ICGC_NaiveOnly_Samples.csv", row.names = FALSE)

############################################ OUTDATED ##############################################################

##### Identify Chemo Naive Samples OUTDATED

library(dplyr)
library(purrr)
library(readr)

#set working directory
setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")

#Identify all files
file_path <- "/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Raw_Clinical_Data/all_projects"
tsv_files <- list.files(file_path, pattern = "\\.tsv$", full.names = TRUE)

#Extract all samples with treatment information
combined_data <- tsv_files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE, col_types = cols(.default = "c")), .id = "file_name")

combined_data <- combined_data[,-1]

# Add colnames  
header <- read_tsv("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only/ICGC_Raw_Clinical_Data/donor_therapy_header.tsv")
colnames(combined_data) <- colnames(header)

# Extract Chemo Naive 

naive_samples <- combined_data[combined_data$first_therapy_type == "no treatment",] # 1024 samples in total across 
naive_samples <- naive_samples[,colnames(naive_samples) %in% c("icgc_donor_id", "project_code", 
                                                               "submitted_donor_id", "first_therapy_type")]

# Save Samples used

write.csv(naive_samples, file = "ICGC_NaiveOnly_Samples.csv", row.names = FALSE)

