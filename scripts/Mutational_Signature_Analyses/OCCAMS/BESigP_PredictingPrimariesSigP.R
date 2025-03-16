# can BE Sig17+ predict whether the matched primary is BELike.Sig17+?

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier")
library(dplyr)
library(epitools)

# Load Pairs from Maria
load("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/BarrettsNaiveComparison/pairs.RData")

# Extract the OCCAMS_ID for only matched Primary Tumors and Barretts from the paired annotations
pairs <- pairs %>%
  group_by(OCCAMS_ID) %>%
  filter(any(Category == 'PrimaryTumour') & any(Category == "Barretts")) %>%
  ungroup() %>%
  distinct(OCCAMS_ID) %>%
  select(OCCAMS_ID)

# Load and extract annotation information about matched samples
load("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/BarrettsNaiveComparison/annotation.sampleIDs.RData")
annotation.sampleIDs <- annotation.sampleIDs[annotation.sampleIDs$OCCAMS_ID %in% pairs$OCCAMS_ID,c("TumourID", "Sample","Category","OCCAMS_ID")]
table(annotation.sampleIDs$OCCAMS_ID) #check
annotation.sampleIDs <- subset(annotation.sampleIDs, select = -c(Sample))

# Load sig17 status annotations (primaries)
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))
names(occams_results_df) <- c("TumourID","Phenotype_Assigned")
occams_results_df <- occams_results_df[occams_results_df$Phenotype_Assigned != "Sig17-",]
occams_results_df$Primary <- ifelse(occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+", "BarrettsLike.Sig17+", "NonBarrettsLike.Sig17+")

# Load sig17 status annotations (barretts)
load("/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata")
ann$TumourID <- rownames(ann)
names(ann) <- c("Phenotype_Assigned","TumourID")
ann$Barretts <- ifelse(ann$Phenotype == "Sig17-", "BE.Sig17-", "BE.Sig17+")


occams_results_df <- subset(occams_results_df, select = -c(Phenotype_Assigned))
ann <- subset(ann, select = -c(Phenotype_Assigned))
occams_results_df <- merge(occams_results_df, ann,
                           by = "TumourID", all = TRUE)

# Collate
matched <- merge(
  occams_results_df,annotation.sampleIDs,
  by = "TumourID"
)

matched <- matched[!matched$OCCAMS_ID %in% c("OCCAMS/PL/056", "OCCAMS/CO/042", "OCCAMS/RS/225"),] #44 samples matched
matched <- matched %>%
  arrange(OCCAMS_ID)

primary_df<- matched[, c("Primary", "OCCAMS_ID")] ; primary_df <- primary_df[!is.na(primary_df$Primary),]
barretts_df <- matched[, c("Barretts", "OCCAMS_ID")] ; barretts_df <- barretts_df[!is.na(barretts_df$Barretts),]
merged_df <- merge(primary_df, barretts_df, by = "OCCAMS_ID")
merged_df$Primary <- factor(merged_df$Primary, levels = c("NonBarrettsLike.Sig17+","BarrettsLike.Sig17+"))


contingency_table <- table(merged_df$Barretts,merged_df$Primary)
epitools::oddsratio.fisher(contingency_table)



# ------ Outdated


# Load Pairs from Maria
load("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/BarrettsNaiveComparison/pairs.RData")

# Extract the OCCAMS_ID for only matched Primary Tumors and Barretts from the paired annotations
pairs <- pairs %>%
  group_by(OCCAMS_ID) %>%
  filter(any(Category == 'PrimaryTumour') & any(Category == "Barretts")) %>%
  ungroup() %>%
  distinct(OCCAMS_ID) %>%
  select(OCCAMS_ID)

# Load and extract annotation information about matched samples
load("/Users/jao/Desktop/MSc_Project/Classifiers_Comparison/BarrettsNaiveComparison/annotation.sampleIDs.RData")
annotation.sampleIDs <- annotation.sampleIDs[annotation.sampleIDs$OCCAMS_ID %in% pairs$OCCAMS_ID,c("TumourID", "Sample","Category","OCCAMS_ID")]
table(annotation.sampleIDs$OCCAMS_ID) #check
annotation.sampleIDs <- subset(annotation.sampleIDs, select = -c(Sample))

# Load sig17 status annotations (primaries)
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))
names(occams_results_df) <- c("TumourID","Phenotype_Assigned")
occams_results_df$Primary <- ifelse(occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+", "BarrettsLike.Sig17+", "NonBarrettsLike")

# Load sig17 status annotations (barretts)
load("/Users/jao/Desktop/MSc_Project/Barretts_Classifier/Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata")
ann$TumourID <- rownames(ann)
names(ann) <- c("Phenotype_Assigned","TumourID")
ann$Barretts <- ifelse(ann$Phenotype == "Sig17-", "BE.Sig17-", "BE.Sig17+")


occams_results_df <- subset(occams_results_df, select = -c(Phenotype_Assigned))
ann <- subset(ann, select = -c(Phenotype_Assigned))
occams_results_df <- merge(occams_results_df, ann,
                           by = "TumourID", all = TRUE)

# Collate

matched <- merge(
  occams_results_df,annotation.sampleIDs,
  by = "TumourID"
)

matched <- matched[!matched$OCCAMS_ID %in% c("OCCAMS/PL/056", "OCCAMS/CO/042", "OCCAMS/RS/225"),] #44 samples matched
matched <- matched %>%
  arrange(OCCAMS_ID)

primary_df<- matched[, c("Primary", "OCCAMS_ID")] ; primary_df <- primary_df[!is.na(primary_df$Primary),]
barretts_df <- matched[, c("Barretts", "OCCAMS_ID")] ; barretts_df <- barretts_df[!is.na(barretts_df$Barretts),]
merged_df <- merge(primary_df, barretts_df, by = "OCCAMS_ID")
merged_df$Primary <- factor(merged_df$Primary, levels = c("NonBarrettsLike","BarrettsLike.Sig17+"))


contingency_table <- table(merged_df$Barretts,merged_df$Primary)
epitools::oddsratio.fisher(contingency_table)

