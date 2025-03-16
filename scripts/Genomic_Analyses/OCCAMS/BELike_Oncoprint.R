# ----- Oncoprint for OCCAMS cohort EAC drivers

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic Analysis/genomic_changes/Oncoprint")
library(GenomicRanges)
library(ComplexHeatmap)
library(reshape2)

drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

# ----- BarrettsLike.Sig17+ Samples

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
BELikeSamples <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"]

# - Generate CNV matrix

# Load OAC CNV data
oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$SampleID <- sub("_vs_.*", "", oac_cna$sample)

# Extract only BELike samples
BELike_cna <- oac_cna[which(oac_cna$SampleID %in% BELikeSamples),]

# Load OAC driver gene locations
hg38 <- read.delim("hg38_geneLocations.txt", header = TRUE, sep = "\t")

gr.hg38 <- GRanges(
  seqnames = hg38$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38$Gene.start..bp., end=hg38$Gene.end..bp.),
  gene=hg38$Gene.name)

gr.BElike.cns <- GRanges(
  seqnames = BELike_cna$chromosome,
  ranges=IRanges(start=BELike_cna$start, end=BELike_cna$end),
  CN=BELike_cna$CNchange,
  sample=BELike_cna$SampleID)

# Find Gene corresponding to the CNA

BELike.ovlp <- data.frame(findOverlaps(gr.BElike.cns,gr.hg38, type="any"))

gr.BELike.overlap <- cbind(as.data.frame(gr.BElike.cns[as.matrix(BELike.ovlp)[,1]]),
                      as.data.frame(gr.hg38[as.matrix(BELike.ovlp)[,2]])) 

gr.BELike.overlap$AMP <- sapply(gr.BELike.overlap$CN, function(x) ifelse(x=="AMP",1,0))
gr.BELike.overlap$DEL <- sapply(gr.BELike.overlap$CN, function(x) ifelse(x=="DEL",1,0))

mat.BELike.cn.amp <- acast(gene~sample, value.var = "AMP",fun.aggregate = sum,
                      data=gr.BELike.overlap)
mat.BELike.cn.amp[mat.BELike.cn.amp > 1] <- 1
mat.BELike.cn.del <- acast(gene~sample, value.var = "DEL",fun.aggregate = sum,
                      data=gr.BELike.overlap)
mat.BELike.cn.del[mat.BELike.cn.del > 1] <- 1

# - Generate SNV and Indel matrix

# Load OAC SNV data 
oac_SNVs.nonsynonymous <- read.delim("SNVs.nonsynonymous.OAC.txt")
oac_SNVs.nonsynonymous$SampleID <- sub("_vs_.*", "", oac_SNVs.nonsynonymous$Sample)

# Load OAC Indel data
oac_indels <- read.delim("indels.OAC.txt")
oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample)  

# Extract only BELike.Sig17 Samples
BELike_snv <- oac_SNVs.nonsynonymous[which(oac_SNVs.nonsynonymous$SampleID %in% BELikeSamples),]
BELike_indel <- oac_indels[which(oac_indels$SampleID %in% BELikeSamples),]

snvs.BELike.selected <- unique(BELike_snv[which(BELike_snv$Gene %in% drivers),])
indels.BELike.selected <- unique(BELike_indel[which(BELike_indel$Gene %in% drivers),])

mat.BELike.snvs <- array(0,c(length(drivers),length(unique(c(snvs.BELike.selected$SampleID,indels.BELike.selected$SampleID,colnames(mat.BELike.cn.amp))))))
rownames(mat.BELike.snvs) <- drivers
colnames(mat.BELike.snvs) <- unique(c(snvs.BELike.selected$SampleID,indels.BELike.selected$SampleID,colnames(mat.BELike.cn.amp)))
for (i in 1:nrow(snvs.BELike.selected)) {
  print(i)
  mat.BELike.snvs[snvs.BELike.selected[i,]$Gene,snvs.BELike.selected[i,]$SampleID] <- 1
}

mat.BELike.indels <- array(0,c(length(drivers),length(unique(c(snvs.BELike.selected$SampleID,indels.BELike.selected$SampleID,colnames(mat.BELike.cn.amp))))))
rownames(mat.BELike.indels) <- drivers
colnames(mat.BELike.indels) <- unique(c(snvs.BELike.selected$SampleID,indels.BELike.selected$SampleID,colnames(mat.BELike.cn.amp)))
for (i in 1:nrow(indels.BELike.selected)) {
  print(i)
  mat.BELike.indels[indels.BELike.selected[i,]$Gene,indels.BELike.selected[i,]$SampleID] <- 1
}


# - Generate Mat_list input for oncoprint heatmap

mat_list.BELike <- list()
mat_list.BELike$snv <- mat.BELike.snvs[drivers,BELikeSamples]
mat_list.BELike$indel <- mat.BELike.indels[drivers,BELikeSamples]
mat.BELike.cn.amp[drivers, BELikeSamples][is.na(mat.BELike.cn.amp[drivers, BELikeSamples])] <- 0
mat_list.BELike$amp <- mat.BELike.cn.amp[drivers,BELikeSamples]
mat.BELike.cn.del[drivers, BELikeSamples][is.na(mat.BELike.cn.del[drivers, BELikeSamples])] <- 0
mat_list.BELike$del <- mat.BELike.cn.del[drivers,BELikeSamples]
names(mat_list.BELike) <- c("SNV","Indel","AMP","DEL")

save(mat_list.BELike, file = "mat_list.BELike.Rdata")

# - Generate Oncoprint Heatmap

# Load Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0
mutsig.BELike <- primaries_sigs_complete[BELikeSamples,]

altered_nums<- Reduce( "+" , lapply(mat_list.BELike, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.BELike<-  lapply(mat_list.BELike, function(x) x[slice_indx,])


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


heatmap_legend_param = list(title = "Alterations", at = c("SNV","Indel","AMP","DEL"),
                            labels = c("Non-synonymous Mutation","Indel", "Amplification", "Deletion"))

pdf("oncoprint.drivers.BELike.pdf",h=6, w=11)
oncoPrint(mat_list.BELike,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS BeLike.Sig17+ Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.BELike[,"SBS17a"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"))),  # Adjusted label size
            SBS17b = anno_barplot(mutsig.BELike[,"SBS17b"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS2 = anno_barplot(mutsig.BELike[,"SBS2"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS8 = anno_barplot(mutsig.BELike[,"SBS8"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS18 = anno_barplot(mutsig.BELike[,"SBS18"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS35 = anno_barplot(mutsig.BELike[,"SBS35"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS40 = anno_barplot(mutsig.BELike[,"SBS40"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS41 = anno_barplot(mutsig.BELike[,"SBS41"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS44 = anno_barplot(mutsig.BELike[,"SBS44"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          ))
dev.off()
