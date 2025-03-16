setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Survival_Analysis")
library(dplyr)
library(survival)
library(forestplot)

# ----- Load and format inputs

# Load Filtered and TMM and CPM Normalized Inputs and annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/Processed_ExpressionCounts.Rdata")
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/Processed_ExpressionCounts_Annotations.Rdata")

# Only use robust differentially expressed gene list
dea <- read.delim("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/Robust_commondegs.txt", header = F)$V1
#dea <- dea[!dea %in% c("GALNT9", "CYP1A1","PNMA5")]
count_norm <- count_norm[rownames(count_norm) %in% dea,]

# Only use Naivelike.Sig17+ and Sig17- samples
annotations <- data.frame(
  subgroup = annotations,
  row.names = names(annotations))
annotations <- annotations[annotations$subgroup %in% c("NaiveLike.Sig17+", "BarrettsLike.Sig17+", "TreatedLike.Sig17+"), , drop = FALSE]
annotations$subgroup <- as.character(annotations$subgroup)

count_norm <- as.data.frame(t(count_norm[,colnames(count_norm) %in% rownames(annotations)]))
all(rownames(annotations) == rownames(count_norm)) #check
count_norm$TumourID <- rownames(count_norm)

# Load and recode clinical data
clin <- read.csv("clin_data_20220609.csv")

# Read in stratified samples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))
occams_results_df <- occams_results_df[occams_results_df$Phenotype_Assigned %in% c("BarrettsLike.Sig17+", "NaiveLike.Sig17+"),]
#occams_results_df$Phenotype_Assigned <- ifelse(occams_results_df$Phenotype_Assigned == "Sig17-", "Sig17-", "Sig17+")
#occams_results_df <- occams_results_df[occams_results_df$Phenotype_Assigned == "Sig17+",]

clin <- clin %>% #279 samples (156 Sig17+) (133 BElike and Naivelike)
  filter(clin$TumourID %in% rownames(count_norm) & clin$TumourID %in% occams_results_df$Sample) %>%
  mutate(status = recode(Patient.Died.c, 'no' = 0, 'yes'=1)) %>%
  select(TumourID, Weeks.Survival.c, status) %>%
  distinct(TumourID, .keep_all = TRUE) 

clin <- inner_join(clin, count_norm, by = "TumourID")
rownames(clin) <- clin$TumourID ; clin <- subset(clin, select = -c(TumourID))

# Perform Cox regression for each gene
cox_results <- lapply(dea, function(gene) {
  formula <- as.formula(paste("Surv(Weeks.Survival.c, status) ~", gene))
  cox_model <- coxph(formula, data = clin)
  
  # Return summary information as a data frame
  data.frame(
    Gene = gene,
    HR = exp(coef(cox_model)[gene]),  # Hazard ratio for the gene
    Lower.CI = exp(confint(cox_model)[gene, 1]),  # Lower confidence interval
    Upper.CI = exp(confint(cox_model)[gene, 2]),  # Upper confidence interval
    P.value = summary(cox_model)$coefficients[gene, 5]  # P-value for the gene
  )
})


# Combine results into a single data frame
results_table <- do.call(rbind, cox_results)
results_table$padj <- p.adjust(results_table$P.value, method = "BH")
results_table <- results_table %>% arrange(padj)

# ----- Forest Plot
# Prepare labeltext (including the header row)
labeltext <- cbind(
  c("Gene", results_table$Gene),  
  c("q-value", formatC(results_table$padj, format = "f", digits = 3)), 
  c("HR (95% CI)", paste0(round(results_table$HR, 2), " (", 
                          round(results_table$Lower.CI, 2), "-", 
                          round(results_table$Upper.CI, 2), ")")))

# Generate forest plot
pdf("occams_cox_forestplot.pdf", width = 6, height = 8)
forestplot(
  labeltext = labeltext,
  mean = c(NA, results_table$HR),         
  lower = c(NA, results_table$Lower.CI),    
  upper = c(NA, results_table$Upper.CI),    
  zero = 1,                               
  xlab = "Hazard Ratio (HR)",
  graph.pos = "left",
  txt_gp = fpTxtGp(label = gpar(fontsize = 10), xlab = gpar(fontsize = 15), ticks = gpar(fontsize = 15)),
  xticks.digits = 1)
dev.off()

