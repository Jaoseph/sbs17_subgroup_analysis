##### Identify Samples with no Indel Deletions
setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")

# Load Library
library(deconstructSigs)

# Load hg19 signatures (for ICGC samples)
signatures.sbs96_hg19 <- read.table('MutationalSignatureReferences/COSMIC_v3.3.1_SBS_GRCh37.txt', h=T)
sig_sbs96_types <- signatures.sbs96_hg19$Type
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, 2:ncol(signatures.sbs96_hg19)]
signatures.sbs96_hg19 <- data.frame(t(signatures.sbs96_hg19))
colnames(signatures.sbs96_hg19) <- sig_sbs96_types
signatures.sbs96_hg19 <- signatures.sbs96_hg19[, colnames(signatures.cosmic)]

# Load COSMIC indel signatures
signatures.id83 <- read.table('MutationalSignatureReferences/COSMIC_v3.3_ID_GRCh37.txt', h=T)
sig_id83_types <- signatures.id83$Type
signatures.id83 <- signatures.id83[, 2:ncol(signatures.id83)]
signatures.id83 <- data.frame(t(signatures.id83))
colnames(signatures.id83) <- sig_id83_types

# Load in datasets and tweak into deconstructSigs inputs
load('ICGC_Naive_mt_tally.Rdata')

sigs.sbs96_input <- as.data.frame(mt_tally.icgc_naive$SBS_96)
sigs.sbs96_input <- sigs.sbs96_input[,colnames(signatures.sbs96_hg19)]

sigs.id83_input <- as.data.frame(mt_tally.icgc_naive$ID_83)
sigs.id83_input <- sigs.id83_input[,colnames(signatures.id83)]

#Check if there is NA in sigs_input for deconstructSigs

#Check NA values
if(any(is.na(sigs.id83_input))) {
  print("NA values found in ID83 input")
} else {
  print("No NA values found in ID83 input")
}

#Check if there is any rows that do not sum to 1

sigs.id83_input_normalized <- as.data.frame(t(apply(sigs.id83_input, 1, function(x) {
  total <- sum(x)
  if(total == 0) {
    # If the total is 0, return a row of zeros to avoid division by zero
    return(rep(0, length(x)))
  } else {
    # Otherwise, return the proportion of each count relative to the row total
    return(x / total)
  }
})))

row_sums_normalized <- rowSums(sigs.id83_input_normalized)
if(all(abs(row_sums_normalized - 1) < .Machine$double.eps*100)) {
  print("Each row now sums to 1.")
} else {
  print("Some rows do not sum to 1 after normalization.")
}

# Ensure that sigs.id83_input is a data frame
sigs.id83_input <- data.frame(sigs.id83_input)

# Check that it is now a data frame
is.data.frame(sigs.id83_input)

rows_not_summing_to_one <- which(abs(row_sums_normalized - 1) > .Machine$double.eps*100)
#print("Rows not summing to 1 after normalization:")
#print(rows_not_summing_to_one)

if (length(rows_not_summing_to_one) > 0) {
  for (i in rows_not_summing_to_one) {
    cat("Details for row", i, ":\n")
    # Convert row to a character vector and collapse it into a single string to print
    original_counts <- paste(sigs.id83_input[i, ], collapse = ", ")
    normalized_counts <- paste(sigs.id83_input_normalized[i, ], collapse = ", ")
    cat("Original counts: ", original_counts, "\n")
    cat("Normalized proportions: ", normalized_counts, "\n")
    cat("Sum of normalized proportions: ", row_sums_normalized[i], "\n\n")
  }
} else {
  cat("All rows sum to 1 after normalization.\n")
}

# Additionally, check for any NA or NaN values after normalization
any_nans <- any(sapply(sigs.id83_input_normalized, function(x) { any(is.nan(x)) }))
cat("Any NaN values after normalization: ", any_nans, "\n")

row_names_not_summing_to_one <- rownames(sigs.id83_input_normalized)[rows_not_summing_to_one]

# Print the vector of row names not summing to one
print(row_names_not_summing_to_one)

rows_to_extract <- sigs.id83_input[rows_not_summing_to_one, ]
print(rows_to_extract)
save(rows_to_extract, file = "RowsNotSummingToOneDetails.RData")

no_ID_samples <- rownames(rows_to_extract)

write.table(no_ID_samples, file = 'no_ID_samples.txt')
no_ID_samples <- read.table('no_ID_samples.txt', h=T)
no_ID_samples <- no_ID_samples$Samples
