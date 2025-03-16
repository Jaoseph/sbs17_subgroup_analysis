setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
# Define the list of required packages
packages <- c("doSNOW", "e1071", "parallel", "snow", "tibble", "dplyr", "readr")

# Install all packages from CRAN
install.packages(packages, repos = "https://cran.r-project.org")

library(parallel)
library(snow)
library(doSNOW)
library(tibble)
library(dplyr)
library(e1071)
library(readr)

print("Starting")
# ---- 2. Model Selection
source("R/train_predict_functions.R")

# Load data annotations + Feature mapped
load("sysSVM2_input.Rdata")

# Separate the training and prediction sets (i.e. canonical drivers, and the rest of genes), and perform data normalisation
canonical_drivers = readRDS("example_data/canonical_drivers.rds")
sysSVM_data = prepare_trainingPrediction(sysSVM2_input, canonical_drivers, output_dir = "/mnt/data/jao/sysSVM/sysSVM2/test_sysSVM2")

cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Running cv_stats\n")
cv_stats = run_crossValidation_par(iters = 1000,
                                   cores = 70, 
                                   inPath = "~/test_sysSVM2",
                                   outPath = "~/test_sysSVM2",
                                   parallelLib = "parallel")

saveRDS(cv_stats, file = "~/test_sysSVM2/cv_stats.rds")

# Identify the best parameter combinations
print("model selection")
model_selection = selectParams_from_CVstats(cv_stats, output_dir = "~/test_sysSVM2")

saveRDS(model_selection, file = "~/test_sysSVM2/model_selection.rds")
