# ----- Oncoprint for OCCAMS cohort EAC drivers

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic Analysis/genomic_changes/Oncoprint")
library(GenomicRanges)
library(ComplexHeatmap)
library(reshape2)

drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

# ----- BarrettsLike.Sig17+ Samples

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
NegativeSamples <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "Sig17-"]

# - Generate CNV matrix

# Load OAC CNV data
oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$SampleID <- sub("_vs_.*", "", oac_cna$sample)

# Extract only Negative samples
Negative_cna <- oac_cna[which(oac_cna$SampleID %in% NegativeSamples),]

# Load OAC driver gene locations
hg38 <- read.delim("hg38_geneLocations.txt", header = TRUE, sep = "\t")

gr.hg38 <- GRanges(
  seqnames = hg38$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38$Gene.start..bp., end=hg38$Gene.end..bp.),
  gene=hg38$Gene.name)

gr.BElike.cns <- GRanges(
  seqnames = Negative_cna$chromosome,
  ranges=IRanges(start=Negative_cna$start, end=Negative_cna$end),
  CN=Negative_cna$CNchange,
  sample=Negative_cna$SampleID)

# Find Gene corresponding to the CNA

Negative.ovlp <- data.frame(findOverlaps(gr.BElike.cns,gr.hg38, type="any"))

gr.Negative.overlap <- cbind(as.data.frame(gr.BElike.cns[as.matrix(Negative.ovlp)[,1]]),
                           as.data.frame(gr.hg38[as.matrix(Negative.ovlp)[,2]])) 

gr.Negative.overlap$AMP <- sapply(gr.Negative.overlap$CN, function(x) ifelse(x=="AMP",1,0))
gr.Negative.overlap$DEL <- sapply(gr.Negative.overlap$CN, function(x) ifelse(x=="DEL",1,0))

mat.Negative.cn.amp <- acast(gene~sample, value.var = "AMP",fun.aggregate = sum,
                           data=gr.Negative.overlap)
mat.Negative.cn.amp[mat.Negative.cn.amp > 1] <- 1
mat.Negative.cn.del <- acast(gene~sample, value.var = "DEL",fun.aggregate = sum,
                           data=gr.Negative.overlap)
mat.Negative.cn.del[mat.Negative.cn.del > 1] <- 1

# - Generate SNV and Indel matrix

# Load OAC SNV data 
oac_SNVs.nonsynonymous <- read.delim("SNVs.nonsynonymous.OAC.txt")
oac_SNVs.nonsynonymous$SampleID <- sub("_vs_.*", "", oac_SNVs.nonsynonymous$Sample)

# Load OAC Indel data
oac_indels <- read.delim("indels.OAC.txt")
oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample)  

# Extract only Negative.Sig17 Samples
Negative_snv <- oac_SNVs.nonsynonymous[which(oac_SNVs.nonsynonymous$SampleID %in% NegativeSamples),]
Negative_indel <- oac_indels[which(oac_indels$SampleID %in% NegativeSamples),]

snvs.Negative.selected <- unique(Negative_snv[which(Negative_snv$Gene %in% drivers),])
indels.Negative.selected <- unique(Negative_indel[which(Negative_indel$Gene %in% drivers),])

mat.Negative.snvs <- array(0,c(length(drivers),length(unique(c(snvs.Negative.selected$SampleID,indels.Negative.selected$SampleID,colnames(mat.Negative.cn.amp))))))
rownames(mat.Negative.snvs) <- drivers
colnames(mat.Negative.snvs) <- unique(c(snvs.Negative.selected$SampleID,indels.Negative.selected$SampleID,colnames(mat.Negative.cn.amp)))
for (i in 1:nrow(snvs.Negative.selected)) {
  print(i)
  mat.Negative.snvs[snvs.Negative.selected[i,]$Gene,snvs.Negative.selected[i,]$SampleID] <- 1
}

mat.Negative.indels <- array(0,c(length(drivers),length(unique(c(snvs.Negative.selected$SampleID,indels.Negative.selected$SampleID,colnames(mat.Negative.cn.amp))))))
rownames(mat.Negative.indels) <- drivers
colnames(mat.Negative.indels) <- unique(c(snvs.Negative.selected$SampleID,indels.Negative.selected$SampleID,colnames(mat.Negative.cn.amp)))
for (i in 1:nrow(indels.Negative.selected)) {
  print(i)
  mat.Negative.indels[indels.Negative.selected[i,]$Gene,indels.Negative.selected[i,]$SampleID] <- 1
}


# - Generate Mat_list input for oncoprint heatmap
NegativeSamples <- NegativeSamples[NegativeSamples %in% colnames(mat.Negative.cn.amp)]

mat_list.Negative <- list()
mat_list.Negative$snv <- mat.Negative.snvs[drivers,NegativeSamples]
mat_list.Negative$indel <- mat.Negative.indels[drivers,NegativeSamples]
mat.Negative.cn.amp[drivers, NegativeSamples][is.na(mat.Negative.cn.amp[drivers, NegativeSamples])] <- 0
mat_list.Negative$amp <- mat.Negative.cn.amp[drivers,NegativeSamples]
mat.Negative.cn.del[drivers, NegativeSamples][is.na(mat.Negative.cn.del[drivers, NegativeSamples])] <- 0
mat_list.Negative$del <- mat.Negative.cn.del[drivers,NegativeSamples]
names(mat_list.Negative) <- c("SNV","Indel","AMP","DEL")

save(mat_list.Negative, file = "mat_list.Negative.Rdata")

# - Generate Oncoprint Heatmap

# Load Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0
mutsig.Negative <- primaries_sigs_complete[NegativeSamples,]

altered_nums<- Reduce( "+" , lapply(mat_list.Negative, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.Negative<-  lapply(mat_list.Negative, function(x) x[slice_indx,])


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
                            labels = c("Non-synonymous Mutation","Indel", "Amplification", "Deletion"))

pdf("oncoprint.drivers.Negative.pdf",h=6,w=11)
oncoPrint(mat_list.Negative,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS Sig17- Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.Negative[,"SBS17a"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"))),  # Adjusted label size
            SBS17b = anno_barplot(mutsig.Negative[,"SBS17b"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS2 = anno_barplot(mutsig.Negative[,"SBS2"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS8 = anno_barplot(mutsig.Negative[,"SBS8"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS18 = anno_barplot(mutsig.Negative[,"SBS18"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS35 = anno_barplot(mutsig.Negative[,"SBS35"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS40 = anno_barplot(mutsig.Negative[,"SBS40"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS41 = anno_barplot(mutsig.Negative[,"SBS41"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS44 = anno_barplot(mutsig.Negative[,"SBS44"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          ))

dev.off()


