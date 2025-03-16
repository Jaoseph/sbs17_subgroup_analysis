setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Survival_Analysis")
library(dplyr)
library(survival)

# Load and recode clinical data
tcga_annotation <- read.delim('esca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep="\t")
tcga_annotation <- tcga_annotation[tcga_annotation$Cancer.Type.Detailed == "Esophageal Adenocarcinoma",]
tcga_annotation <- tcga_annotation %>%
  mutate(status = recode(Overall.Survival.Status, "0:LIVING" = 0, "1:DECEASED" = 1)) %>%
  select(Patient.ID,Overall.Survival..Months.,status)

# Load TCGA mutational Signature 17 group classifications

load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga <- subset(ann_tcga, select = Phenotype_Assigned)
ann_tcga$Patient.ID <- rownames(ann_tcga)
tcga_complete <- merge(ann_tcga, tcga_annotation, by = "Patient.ID")
tcga_complete <- tcga_complete[tcga_complete$Patient.ID != "TCGA-L5-A4OT",]

tcga_complete$Phenotype_Assigned <- factor(tcga_complete$Phenotype_Assigned, 
                                  levels = c("NaiveLike.Sig17+", "Sig17-"))

# Generate Surv Object
surv <- survfit(Surv(Overall.Survival..Months., status) ~ Phenotype_Assigned, data = tcga_complete)

my_colors <- c(
  "NaiveLike.Sig17+" = wes_palette("AsteroidCity3")[4],
  "Sig17-" = wes_palette("AsteroidCity3")[2])

# Plot Survival differences between sub-groups
pdf("tcga_Subtype_SurvDiff.pdf", w=10,h=10)
ggsurvplot(
  surv,
  data = tcga_complete,
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
  legend.labs = c("NaiveLike.Sig17+", "Sig17-"))
dev.off()
