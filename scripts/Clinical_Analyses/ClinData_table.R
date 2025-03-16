setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Clinical_Analysis")

library(dplyr)

# Load TCGA clinical data
tcga_clin <- read.delim('esca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep="\t")
tcga_clin <- tcga_clin[tcga_clin$Cancer.Type.Detailed == "Esophageal Adenocarcinoma",]
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")

# Remove BarrettsLike.Sig17+ sample
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",]
tcga_clin <- tcga_clin[tcga_clin$Patient.ID %in% rownames(ann_tcga),]

# Load OCCCAMS clinical data
occams_clin <- read.csv("/Users/jao/Desktop/MSc_Project/Clinical_Data/clin_data_20220609.csv")
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_clin <- occams_clin[occams_clin$TumourID %in% occams_results_df$Sample,]

# Re-format clinical data
occams_clin$Patient.Died.c <- ifelse(occams_clin$Patient.Died.c == "no", "0:LIVING", "1:DECEASED")
occams_clin$DI.PatientGender <- ifelse(occams_clin$DI.PatientGender == "male", "Male", "Female")

# Standardize column

tcga_clin <- tcga_clin %>%
  rename(
    SampleID = Patient.ID,
    DiagnosisAge = Diagnosis.Age,
    Gender = Sex,
    SurvivalStatus = Overall.Survival.Status,
    TumorStage = American.Joint.Committee.on.Cancer.Tumor.Stage.Code,
    LymphNodeStage = Neoplasm.Disease.Lymph.Node.Stage.American.Joint.Committee.on.Cancer.Code,
    MetastasisStage = American.Joint.Committee.on.Cancer.Metastasis.Stage.Code
  ) %>%
  mutate(
    SurvivalWeeks = Overall.Survival..Months. * 4.33,
    Dataset = "TCGA"
  )

occams_clin <- occams_clin %>%
  rename(
    SampleID = TumourID,
    DiagnosisAge = DI.ageAtDiagnosis,
    Gender = DI.PatientGender,
    SurvivalWeeks = Weeks.Survival.c,
    SurvivalStatus = Patient.Died.c,
    TumorStage = PS.TStage.PrimaryTumour.FinalPretreatmentStaging,
    LymphNodeStage = PS.NStage.PrimaryTumour.FinalPretreatmentStaging.TNM7,
    MetastasisStage = PS.MStage.PrimaryTumour.FinalPretreatmentStaging,
    MandardScore = RP.MandardScoreForResponse
  ) %>%
  mutate(
    Dataset = "OCCAMS"
  )

# Common clinical data
common_clin <- c("SampleID","DiagnosisAge", "Gender", "SurvivalWeeks", "SurvivalStatus", 
                 "TumorStage", "LymphNodeStage", "MetastasisStage", "Dataset")


tcga_common <- tcga_clin %>% select(all_of(common_clin))
occams_common <- occams_clin %>% select(all_of(common_clin))

overall_clinical <- bind_rows(tcga_common, occams_common)

# Re-format clinical data
overall_clinical <- overall_clinical %>%
  mutate(
    # Convert to uppercase to standardize (already done)
    TumorStage = toupper(TumorStage),
    LymphNodeStage = toupper(LymphNodeStage),
    MetastasisStage = toupper(MetastasisStage),
    
    # Merge subcategories into broader categories
    TumorStage = case_when(
      TumorStage %in% c("T1A", "T1B") ~ "T1",
      TumorStage %in% c("T4A", "T4B") ~ "T4",
      TumorStage == "TIS" ~ "T0",  
      TRUE ~ TumorStage  
    ),
    
    MetastasisStage = case_when(
      MetastasisStage == "M1A" ~ "M1",
      TRUE ~ MetastasisStage
    )
  )

# Create table
library(tableone)

# Define variables to summarize
categorical_vars <- c("Gender", "SurvivalStatus", "TumorStage", "LymphNodeStage", "MetastasisStage")
continuous_vars <- c("DiagnosisAge", "SurvivalWeeks")

# Create a summary table comparing TCGA and OCCAMS
summary_table <- CreateTableOne(vars = c(continuous_vars, categorical_vars), 
                                strata = "Dataset",  # Compare by dataset
                                data = overall_clinical, 
                                factorVars = categorical_vars,
                                testExact = TRUE)
summary_df <- print(summary_table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE) %>%
  as.data.frame()

library(openxlsx)
write.xlsx(summary_df, "Cohort_ComparisonSummary_table.xlsx", sheetName = "Table S", rowNames=TRUE)


# ------ Compare by Sig17 subgroup? (test)

load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))
names(occams_results_df) <- c("SampleID", "Phenotype_Assigned")
occams_results_df <- merge(occams_results_df, occams_clin, by = "SampleID")


occams_results_df <- subset(occams_results_df, select = c(SampleID,Phenotype_Assigned,Gender,SurvivalStatus,
                                                          TumorStage,LymphNodeStage,MetastasisStage,DiagnosisAge,
                                                          SurvivalWeeks,MandardScore))

# Re-format clinical data
occams_results_df <- occams_results_df %>%
  mutate(
    # Convert to uppercase to standardize (already done)
    TumorStage = toupper(TumorStage),
    LymphNodeStage = toupper(LymphNodeStage),
    MetastasisStage = toupper(MetastasisStage),
    
    # Merge subcategories into broader categories
    TumorStage = case_when(
      TumorStage %in% c("T1A", "T1B") ~ "T1",
      TumorStage %in% c("T4A", "T4B") ~ "T4",
      TumorStage == "TIS" ~ "T0",  
      TRUE ~ TumorStage  
    ),
    
    MetastasisStage = case_when(
      MetastasisStage == "M1A" ~ "M1",
      TRUE ~ MetastasisStage
    )
  )

# Define variables to summarize
OCCAMS_categorical_vars <- c("Gender", "SurvivalStatus", "TumorStage", "LymphNodeStage", "MetastasisStage", "MandardScore")
OCCAMS_continuous_vars <- c("DiagnosisAge", "SurvivalWeeks")


# Create a summary table comparing all 4 subtypes
summary_table <- CreateTableOne(vars = c(OCCAMS_continuous_vars, OCCAMS_categorical_vars), 
                                strata = "Phenotype_Assigned",  # Compare by Subgroup
                                data = occams_results_df, 
                                factorVars = OCCAMS_categorical_vars)
summary_df <- print(summary_table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE) %>%
  as.data.frame()

library(rcompanion)
pairwise_TStage <- pairwiseNominalIndependence(
  table(occams_results_df$Phenotype_Assigned, occams_results_df$TumorStage),
  fisher = TRUE,
  gtest  = FALSE,
  chisq  = FALSE)
pairwise_TStage <- pairwise_TStage[!(pairwise_TStage$Comparison %in% c("BarrettsLike.Sig17+ : TreatedLike.Sig17+", "Sig17- : TreatedLike.Sig17+")),] 
pairwise_TStage$p.adj.Fisher <- p.adjust(pairwise_TStage$p.Fisher, method = "BH")
pairwise_TStage

#> pairwise_TStage
#Comparison p.Fisher p.adj.Fisher
#1 BarrettsLike.Sig17+ : NaiveLike.Sig17+   0.1360        0.145
#2           BarrettsLike.Sig17+ : Sig17-   0.1450        0.145
#4              NaiveLike.Sig17+ : Sig17-   0.0215        0.086
#5  NaiveLike.Sig17+ : TreatedLike.Sig17+   0.0825        0.145


# ---- Trial

occams_results_df_2 <- subset(occams_results_df, select = c(SampleID,Phenotype_Assigned,Gender,SurvivalStatus,
                                                          TumorStage,LymphNodeStage,MetastasisStage,DiagnosisAge,
                                                          SurvivalWeeks,MandardScore))
occams_results_df_2 <- occams_results_df_2[occams_results_df_2$Phenotype_Assigned %in% c("BarrettsLike.Sig17+", "NaiveLike.Sig17+"),]

# Create a summary table comparing BElike and NaiveLike subtypes
summary_table <- CreateTableOne(vars = c(OCCAMS_continuous_vars, OCCAMS_categorical_vars), 
                                strata = "Phenotype_Assigned",  # Compare by Subgroup
                                data = occams_results_df_2, 
                                factorVars = OCCAMS_categorical_vars)
summary_df <- print(summary_table, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE) %>%
  as.data.frame()

