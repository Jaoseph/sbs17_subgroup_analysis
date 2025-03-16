##### 
## Creation of signature phenotypes in ICGC Naive_only samples
#####

setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")

### Load libraries
library(mclust)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(writexl)

### Load ICGC deconstructSigs data
load("ICGC_Naive_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata")

### Run mixture modelling using mclust (SBS17-focused clusters)

#sigs_complete (rows = samples, columns = signatures)
set.seed(241202)
#sig_17 <- sigs_complete[,c("SBS17a", "SBS17b")]
#sigs.BIC <- mclustBIC(sig_17, G = 1:2,
#                      modelNames = c('EEI','VII','EII'))
sigs.BIC <- mclustBIC(sigs_complete, G = 1:7,
                      modelNames = c('EEI','VII','EII'))
summary(sigs.BIC)

### Apply optimal clustering in sigs.BIC
mod_sigs.BIC <- Mclust(sigs_complete, x = sigs.BIC)
table(mod_sigs.BIC$classification) # classification

#1    2    3    4    5    6    7 
#1365  220  114   78  374  360  155 


### Form annotation document for visualisation

### Load ICGC sample data (for sampleID matching)
load("/Users/jao/Desktop/MSc_Project/Pan-cancer_Analysis/ICGC_PhenotypeClustering/ICGC/samples.icgc.RData")

samples.icgc <- samples.icgc[samples.icgc$icgc_donor_id %in% rownames(sigs_complete),]
rownames(samples.icgc) <- samples.icgc$icgc_donor_id; samples.icgc <- samples.icgc[,-2, drop = FALSE]
samples.icgc$project_code <- gsub("-.*", "", samples.icgc$project_code)

ann <- data.frame(BIC_clust = factor(mod_sigs.BIC$classification,
                                     levels = 1:length(unique(mod_sigs.BIC$classification))))
ann <- merge(x = ann, y = samples.icgc,
             by = 0, all.x = TRUE)
rownames(ann) <- ann$Row.names; ann <- ann[,-1]
ann <- ann[order(ann$BIC_clust), , drop = FALSE]

#pheno <- c('Sig17-','Sig17+')
#ann$Phenotype <- factor(pheno[as.numeric(ann$BIC_clust)],
#                        levels = sort(pheno))
#ann <- ann[ ,-1]


# Order samples by classification
sigs_order <- as.data.frame(t(sigs_complete[rownames(ann), ]))
sigs_order <- sigs_order[, rownames(ann)]

# Set colours

cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_Pheno <- cols(length(unique(ann$BIC_clust)))
names(cols_Pheno) <- unique(ann$BIC_clust)


cols_cancertype <- c(
  BRCA = "#E6194B", EOPC = "#3CB44B", ESAD = "#FFE119", GACA = "#0082C8",
  ORCA = "#F58231", BTCA = "#911EB4", SKCA = "#46F0F0", MELA = "#F032E6",
  MALY = "#D2F53C", PRAD = "#FABEBE", PACA = "#008080", LAML = "#E6BEFF",
  LICA = "#AA6E28", CLLE = "#FFFAC8", LIAD = "#800000", LIRI = "#AAFFC3",
  LINC = "#808000", OV = "#FFD8B1", RECA = "#000080", ESCA = "#808080",
  CMDI = "#FFFFFF"
)


ann_colors = list(
  Phenotype = cols_Pheno,
  project_code =  cols_cancertype
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

# Generate Heatmap
pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10
         ,filename = 'ICGC_Naive_SignaturesHeatmap_v1.pdf'
)


ann$Phenotype <- ifelse(ann$BIC_clust == 7, "Sig17+", "Sig17-")
ann <- ann[ ,-1, drop = FALSE] #remove BIC_clust
ann <- ann[order(ann$Phenotype), , drop = FALSE]
sigs_order <- sigs_order[, rownames(ann)]

#cols_Pheno <- cols(length(unique(ann$Phenotype)))
cols_Pheno <- c("#E41A1C","#999999")
names(cols_Pheno) <- unique(ann$Phenotype)

ann_colors = list(
  Phenotype = cols_Pheno,
  project_code =  cols_cancertype
)


# Generate Heatmap
pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10
         ,filename = 'ICGC_Naive_SignaturesHeatmap.pdf'
)


save(ann, file = 'ICGC_Naive_IDnormalised_PhenotypeAnnotation_clust.Rdata')


##### counts

result <- ann %>%
  group_by(project_code, Phenotype) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Save
library(openxlsx)
write.xlsx(result, "ICGC_Naive_classification_project_counts.xlsx", sheetName = "Table S", rowNames=TRUE)

total_sig17p <- result %>%
  filter(Phenotype == "Sig17+") %>%
  summarise(Total_Count = sum(Count))

# Print the result
print(total_sig17p) #155
