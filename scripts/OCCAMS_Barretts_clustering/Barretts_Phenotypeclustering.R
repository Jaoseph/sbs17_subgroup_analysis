##### 
## Creation of signature phenotypes in Barretts samples
#####

setwd("/Users/jao/Desktop/MSc_Project/Barretts_Classifier")

### Load libraries
library(mclust)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(writexl)

### Load ICGC deconstructSigs data
load("Barretts_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata")

### Run mixture modelling using mclust (SBS17-focused clusters)

#sigs_complete (rows = samples, columns = signatures)
set.seed(241202)
sigs.BIC <- mclustBIC(sigs_complete, G = 1:4,
                      modelNames = c('EEI','VII','EII'))


summary(sigs.BIC)

### Apply optimal clustering in sigs.BIC
mod_sigs.BIC <- Mclust(sigs_complete, x = sigs.BIC)
table(mod_sigs.BIC$classification) # classification

# 1  2  3  4
#38 65 20  4


### Form annotation document for visualisation

ann <- data.frame(BIC_clust = factor(mod_sigs.BIC$classification,
                                     levels = 1:length(unique(mod_sigs.BIC$classification))))
ann <- ann[order(ann$BIC_clust), , drop =FALSE]

# Order samples by classification

sigs_order <- as.data.frame(t(sigs_complete[rownames(ann), ]))
sigs_order <- sigs_order[, rownames(ann)]

# Set colours

cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_Pheno <- cols(length(unique(ann$BIC_clust)))
names(cols_Pheno) <- unique(ann$BIC_clust)

ann_colors = list(Phenotype = cols_Pheno)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

# Generate Heatmap
pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10
#         ,filename = 'Barretts_SignaturesHeatmap_v1.pdf'
)

ann$Phenotype <- ifelse(ann$BIC_clust == 1, "Sig17+", "Sig17-")
ann <- ann[ ,-1, drop = FALSE] #remove BIC_clust
ann <- ann[order(ann$Phenotype), , drop = FALSE]
sigs_order <- sigs_order[, rownames(ann)]

#cols_Pheno <- cols(length(unique(ann$Phenotype)))
cols_Pheno <- c("#E41A1C","#999999")
names(cols_Pheno) <- unique(ann$Phenotype)

ann_colors = list(Phenotype = cols_Pheno)

# Generate Heatmap
pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10
#         ,filename = 'Barretts_SignaturesHeatmap.pdf'
)

save(ann, file = 'Barretts_IDnormalised_PhenotypeAnnotation_clust.Rdata')


##### counts

result <- ann %>%
  group_by(Phenotype) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Save
library(openxlsx)
write.xlsx(result, "ICGC_Barretts_classification_project_counts.xlsx", sheetName = "Table S", rowNames=TRUE)

total_sig17p <- result %>%
  filter(Phenotype == "Sig17+") %>%
  summarise(Total_Count = sum(Count))

# Print the result
print(total_sig17p) #31
