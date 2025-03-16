# ----- Sub-group Survival analysis

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Survival")
library(dplyr)
library(survival)
library(survminer)
library(ggsurvfit)
library(wesanderson)

# Load clinical data
clin <- read.csv("clin_data_20220609.csv")
clin <- clin %>%
  mutate(status = recode(Patient.Died.c, 'no' = 0, 'yes'=1)) %>%
  select(TumourID, Weeks.Survival.c, status)

# Load classification data
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample, Phenotype_Assigned))
names(occams_results_df) <- c("TumourID", "Phenotype_Assigned")
clin <- merge(clin, occams_results_df, by = "TumourID")

clin$Phenotype_Assigned <- factor(clin$Phenotype_Assigned, 
                                  levels = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-"))
# Rough check 
sum(clin$Phenotype_Assigned == "BarrettsLike.Sig17+") #101
sum(clin$Phenotype_Assigned == "NaiveLike.Sig17+") #249
sum(clin$Phenotype_Assigned == "TreatedLike.Sig17+") #63
sum(clin$Phenotype_Assigned == "Sig17-") #323

# Generate Surv Object
surv <- survfit(Surv(Weeks.Survival.c, status) ~ Phenotype_Assigned, data = clin)

# Plot Survival differences between sub-groups
my_colors <- c(
  "Sig17-" = wes_palette("AsteroidCity3")[2],
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "BarrettsLike.Sig17+" = wes_palette("AsteroidCity3")[3],
  "TreatedLike.Sig17+" = wes_palette("AsteroidCity3")[1])

pdf("Subtype_SurvDiff.pdf", w=10,h=10)
ggsurvplot(
  surv,
  data = clin,
  pval = TRUE,               
  pval.size = 5,             
  pval.method = TRUE,
  pval.method.size = 4,
  conf.int = TRUE,           
  risk.table = TRUE,        
  risk.table.col = "strata", 
  palette = my_colors,
  ggtheme = theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ), 
  xlab = "Survival Weeks",   
  ylab = "Overall survival probability",
  legend.title = "Subgroup",    
  legend.labs = c("BarrettsLike.Sig17+", "NaiveLike.Sig17+", "TreatedLike.Sig17+", "Sig17-"))
dev.off()

