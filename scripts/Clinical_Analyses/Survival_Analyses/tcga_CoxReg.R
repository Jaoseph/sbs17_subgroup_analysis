setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Survival_Analysis")
library(dplyr)
library(survival)

# ----- Load and format inputs

# Load Filtered and TMM and CPM Normalized Inputs and annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/tcga_Processed_ExpressionCounts.Rdata")
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/tcga_Processed_ExpressionCounts_Annotations.Rdata")

# Only use robust differentially expressed gene list
dea <- read.delim("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis/Robust_commondegs.txt", header = F)$V1
count_norm <- count_norm[rownames(count_norm) %in% dea,]

# Only use Naivelike.Sig17+ and Sig17- samples
annotations <- data.frame(
  subgroup = annotations,
  row.names = names(annotations))
annotations$subgroup <- as.character(annotations$subgroup)

count_norm <- as.data.frame(t(count_norm[,colnames(count_norm) %in% rownames(annotations)]))
all(rownames(annotations) == rownames(count_norm)) #check
count_norm <- count_norm[rownames(count_norm) %in% rownames(annotations)[annotations$subgroup == "NaiveLike"], ]
count_norm$Patient.ID <- rownames(count_norm)

# Load and recode clinical data
tcga_annotation <- read.delim('esca_tcga_pan_can_atlas_2018_clinical_data.tsv', sep="\t")
tcga_annotation <- tcga_annotation[tcga_annotation$Cancer.Type.Detailed == "Esophageal Adenocarcinoma",]
tcga_annotation <- tcga_annotation %>%
  mutate(status = recode(Overall.Survival.Status, "0:LIVING" = 0, "1:DECEASED" = 1)) %>%
  select(Patient.ID,Overall.Survival..Months.,status)


tcga <- inner_join(tcga_annotation, count_norm, by = "Patient.ID")
rownames(tcga) <- tcga$Patient.ID ; tcga <- subset(tcga, select = -c(Patient.ID))

# Perform Cox regression for each gene 
cox_results <- lapply(dea, function(gene) {
  formula <- as.formula(paste("Surv(Overall.Survival..Months., status) ~", gene))
  cox_model <- coxph(formula, data = tcga)
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
results_table <- results_table %>% arrange(P.value)

# ----- Forest Plot
# Prepare labeltext (including the header row)
labeltext <- cbind(
  c("Gene", results_table$Gene),  
  c("P-value", formatC(results_table$P.value, format = "f", digits = 3)), 
  c("HR (95% CI)", paste0(round(results_table$HR, 2), " (", 
                          round(results_table$Lower.CI, 2), "-", 
                          round(results_table$Upper.CI, 2), ")")))

# Generate forest plot
pdf("tcga_cox_forestplot.pdf")
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

