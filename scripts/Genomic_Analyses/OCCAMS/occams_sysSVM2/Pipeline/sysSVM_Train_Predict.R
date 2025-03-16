setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/sysSVM2")
source("R/train_predict_functions.R")

# ---- 3. Training

# Load data annotations + Feature mapped
load("sysSVM2_input.Rdata")

# Separate the training and prediction sets (i.e. canonical drivers, and the rest of genes), and perform data normalisation
canonical_drivers = readRDS("example_data/canonical_drivers.rds")
sysSVM_data = prepare_trainingPrediction(sysSVM2_input, canonical_drivers, output_dir = "test_sysSVM2/")

# Load the best parameter combinations
model_selection <- readRDS("test_sysSVM2/model_selection.rds")$best_model_final

# Entire training set is used to train the final sysSVM2 model
trained_sysSVM2 = train_sysSVM2(model_parameters = model_selection, 
                                training_set = sysSVM_data$training_set, 
                                scaling_factors = sysSVM_data$scaling_factors,
                                output_dir = "/test_sysSVM2")


trained_sysSVM <- readRDS("test_sysSVM2/trained_sysSVM.rds")

# ---- 4. Prediction

# trained model can now be used to make predictions on the same cohort
predictions = predict_sysSVM2(trained_sysSVM, 
                              prediction_set = sysSVM_data$prediction_set, 
                              prediction_set_ns = sysSVM_data$prediction_set_ns, 
                              output_dir = "test_sysSVM2")

predictions = readRDS("test_sysSVM2/scores.rds")

# Top-up procedure
drivers_toppedUp = topUp_drivers(all_genes = sysSVM_data,
                                 gene_scores = predictions,
                                 canonical_drivers = canonical_drivers,
                                 n_drivers_per_sample = 10,
                                 output_dir = "test_sysSVM2",
                                 sample_gene_sep = "__")

drivers <- unique(drivers_toppedUp$symbol)
write.table(drivers, file = "sysSVM_drivers.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
