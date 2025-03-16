# ----- Oncoprint for TCGA cohort EAC drivers
setwd('/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Genomic_Analysis')
library(dplyr)
library(ComplexHeatmap)

# ---- Load required annotations and data
# Load Frankel et al. drivers
drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1
drivers[drivers == "H3C2"] <- "HIST1H3B" #H3C2 = HIST1H3B
drivers[drivers == "INSYN2B"] <- "FAM196B" #INSYN2B = FAM196B

# TCGA Sig17 annotations
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/TCGA_Combinedclassification_annotation.Rdata")
ann_tcga <- ann_tcga[ann_tcga$Phenotype_Assigned != "BarrettsLike.Sig17+",] #86 samples

# TCGA mutational data
load("tcga_mutations.Rdata")

# CNV data, extract only EAC samples
cna <- read.delim("esca_tcga_pan_can_atlas_2018/data_cna.txt")
cna <- cna[!duplicated(cna$Hugo_Symbol), ]
rownames(cna) <- cna[,1]
cna <- as.data.frame(t(cna[,-1:-2]))

cna$Patient_ID <- rownames(cna)
cna$Patient_ID <- gsub("\\.", "-" , cna$Patient_ID)
cna$Patient_ID <- sapply(cna$Patient_ID, function(x) substr(x,1,12))
rownames(cna) <- cna$Patient_ID
cna <- cna[which(cna$Patient_ID %in% rownames(ann_tcga)),]
cna <- cna[, -ncol(cna)] ; cna <- as.data.frame(t(cna))

# ---- Transform into matrix

# Amplification and Deletion
cna$gene_name <- rownames(cna)
cna <- cna[cna$gene_name %in% drivers,] ; cna <- cna[,-ncol(cna)] # extract only driver genes

mat.cn.amp <- as.data.frame(ifelse(cna > 1, 1, 0))
mat.cn.del <- as.data.frame(ifelse(cna < -1, 1, 0))

mat.cn.amp <- mat.cn.amp[drivers, , drop = FALSE] # re-order
mat.cn.del <- mat.cn.del[drivers, , drop = FALSE]

# SNV and INDEL
# Extract only driver SNV data (remove silent mutations and select only EAC drivers)
snv <- tcga_mutations[tcga_mutations$Variant_Type == "SNP", c("Hugo_Symbol", "Variant_Classification", 
                                                              "Tumor_Seq_Allele1", "Tumor_Seq_Allele2","Patient_ID")]
snv <- snv[snv$Variant_Classification != "Silent",]
snvs.selected <- unique(snv[which(snv$Hugo_Symbol %in% drivers),])

# Extract only driver INDEL data, (select only EAC drivers)
indels <- tcga_mutations[which(tcga_mutations$Variant_Type %in% c("INS", "DEL")), c("Hugo_Symbol", "Variant_Classification", "Tumor_Seq_Allele1",
                                                                                    "Tumor_Seq_Allele2", "Patient_ID")]
indels.selected <- unique(indels[which(indels$Hugo_Symbol %in% drivers),])

mat.snvs <- array(0,c(length(drivers),length(unique(c(snvs.selected$Patient_ID,indels.selected$Patient_ID,colnames(mat.cn.amp))))))
rownames(mat.snvs) <- drivers
colnames(mat.snvs) <- unique(c(snvs.selected$Patient_ID,indels.selected$Patient_ID,colnames(mat.cn.amp)))
for (i in 1:nrow(snvs.selected)) {
  print(i)
  mat.snvs[snvs.selected[i,]$Hugo_Symbol,snvs.selected[i,]$Patient_ID] <- 1
}

mat.indels <- array(0,c(length(drivers),length(unique(c(snvs.selected$Patient_ID,indels.selected$Patient_ID,colnames(mat.cn.amp))))))
rownames(mat.indels) <- drivers
colnames(mat.indels) <- unique(c(snvs.selected$Patient_ID,indels.selected$Patient_ID,colnames(mat.cn.amp)))
for (i in 1:nrow(indels.selected)) {
  print(i)
  mat.indels[indels.selected[i,]$Hugo_Symbol,indels.selected[i,]$Patient_ID] <- 1
}
# ---- Generate Matrix List for Oncoprint (input)
# NaiveLike.Sig17+ Samples
mat_list.p <- list()
mat_list.p$snv <- mat.snvs[,colnames(mat.snvs) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "NaiveLike.Sig17+"]]
mat_list.p$indel <- mat.indels[,colnames(mat.indels) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "NaiveLike.Sig17+"]]
mat_list.p$amp <- as.matrix(mat.cn.amp[,colnames(mat.cn.amp) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "NaiveLike.Sig17+"]])
mat_list.p$del <- as.matrix(mat.cn.del[,colnames(mat.cn.del) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "NaiveLike.Sig17+"]])
names(mat_list.p) <- c("SNV","Indel","AMP","DEL")

save(mat_list.p, file = "tcga_mat_list.NaiveLike.Rdata")

# Sig17- Samples
mat_list.n <- list()
mat_list.n$snv <- mat.snvs[,colnames(mat.snvs) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "Sig17-"]]
mat_list.n$indel <- mat.indels[,colnames(mat.indels) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "Sig17-"]]
mat_list.n$amp <- as.matrix(mat.cn.amp[,colnames(mat.cn.amp) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "Sig17-"]])
mat_list.n$del <- as.matrix(mat.cn.del[,colnames(mat.cn.del) %in% rownames(ann_tcga)[ann_tcga$Phenotype_Assigned == "Sig17-"]])
names(mat_list.n) <- c("SNV","Indel","AMP","DEL")

save(mat_list.n, file = "tcga_mat_list.Negative.Rdata")

# ---- Generate Oncoprint
# Load Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/TCGA_deconstructSigs/TCGA_deconstructSigs.Rdata")
tcga_sigs_complete[tcga_sigs_complete < 0.05] <- 0

mutsig.NaiveLike <- tcga_sigs_complete[colnames(mat_list.p$SNV),]
mutsig.Negative <- tcga_sigs_complete[colnames(mat_list.n$SNV),]

altered_nums<- Reduce( "+" , lapply(mat_list.p, rowSums))
NUMBER_GENES<- 20
slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]
mat_list.p<-  lapply(mat_list.p, function(x) x[slice_indx,])

altered_nums<- Reduce( "+" , lapply(mat_list.n, rowSums))
NUMBER_GENES<- 20
slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]
mat_list.n<-  lapply(mat_list.n, function(x) x[slice_indx,])


alter_fun_list = list(
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#445E93", col = NA))
  },
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#7EB2DD", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#F93943", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#FFB563", col = NA))
  }
)

col = c(SNV = "#445E93", Indel = "#7EB2DD",
        AMP="#F93943", DEL="#FFB563")


alter_fun_list = list(
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#445E93", col = NA))
  },
  Indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#7EB2DD", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#F93943", col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.9, "mm"), h-unit(1, "mm"), gp = gpar(fill = "#FFB563", col = NA))
  }
)

col = c(SNV = "#445E93", Indel = "#7EB2DD",
        AMP="#F93943", DEL="#FFB563")


#heatmap_legend_param = list(title = "Alterations", at = c("SNV","Indel","AMP","DEL"),
#                            labels = c("Non-synonymous Mutation","Indel", "Amplification", "Deletion"))

pdf("oncoprint.drivers.NaiveLikeTCGA.pdf",h=6, w=11)
oncoPrint(mat_list.p,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS BeLike.Sig17+ Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.NaiveLike[,"SBS17a"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"), at = c(0.1, 0.2))),  # Removed 0
            SBS17b = anno_barplot(mutsig.NaiveLike[,"SBS17b"], ylim = c(0, 0.3), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.15, 0.3))),  # Removed 0
            SBS1 = anno_barplot(mutsig.NaiveLike[,"SBS1"], ylim = c(0, 0.8), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.4, 0.8))),  # Removed 0
            SBS5 = anno_barplot(mutsig.NaiveLike[,"SBS5"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.3, 0.6))),  # Removed 0
            SBS18 = anno_barplot(mutsig.NaiveLike[,"SBS18"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.1, 0.2))),  # Removed 0
            SBS35 = anno_barplot(mutsig.NaiveLike[,"SBS35"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.1, 0.2))),  # Removed 0
            SBS40 = anno_barplot(mutsig.NaiveLike[,"SBS40"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.1, 0.2))),  # Removed 0
            SBS41 = anno_barplot(mutsig.NaiveLike[,"SBS41"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.1, 0.2))),  # Removed 0
            SBS44 = anno_barplot(mutsig.NaiveLike[,"SBS44"], ylim = c(0, 0.2), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 6), at = c(0.1, 0.2))),  # Removed 0
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          )
          
          )
dev.off()

pdf("oncoprint.drivers.NegativeTCGA.pdf",h=6, w=11)
oncoPrint(mat_list.n,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS BeLike.Sig17+ Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.Negative[,"SBS17a"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"))),  # Adjusted label size
            SBS17b = anno_barplot(mutsig.Negative[,"SBS17b"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS1 = anno_barplot(mutsig.Negative[,"SBS1"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS5 = anno_barplot(mutsig.Negative[,"SBS5"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS18 = anno_barplot(mutsig.Negative[,"SBS18"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS35 = anno_barplot(mutsig.Negative[,"SBS35"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS40 = anno_barplot(mutsig.Negative[,"SBS40"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS41 = anno_barplot(mutsig.Negative[,"SBS41"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS44 = anno_barplot(mutsig.Negative[,"SBS44"], ylim = c(0, 0.6), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          ))
dev.off()

