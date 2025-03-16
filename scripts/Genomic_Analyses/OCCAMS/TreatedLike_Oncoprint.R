# ----- Oncoprint for OCCAMS cohort EAC drivers

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic Analysis/genomic_changes/Oncoprint")
library(GenomicRanges)
library(ComplexHeatmap)
library(reshape2)

drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

# ----- BarrettsLike.Sig17+ Samples

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
TreatedLikeSamples <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"]

# - Generate CNV matrix

# Load OAC CNV data
oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$SampleID <- sub("_vs_.*", "", oac_cna$sample)

# Extract only TreatedLike samples
TreatedLike_cna <- oac_cna[which(oac_cna$SampleID %in% TreatedLikeSamples),]

# Load OAC driver gene locations
hg38 <- read.delim("hg38_geneLocations.txt", header = TRUE, sep = "\t")

gr.hg38 <- GRanges(
  seqnames = hg38$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38$Gene.start..bp., end=hg38$Gene.end..bp.),
  gene=hg38$Gene.name)

gr.Treatedlike.cns <- GRanges(
  seqnames = TreatedLike_cna$chromosome,
  ranges=IRanges(start=TreatedLike_cna$start, end=TreatedLike_cna$end),
  CN=TreatedLike_cna$CNchange,
  sample=TreatedLike_cna$SampleID)

# Find Gene corresponding to the CNA

TreatedLike.ovlp <- data.frame(findOverlaps(gr.Treatedlike.cns,gr.hg38, type="any"))

gr.TreatedLike.overlap <- cbind(as.data.frame(gr.Treatedlike.cns[as.matrix(TreatedLike.ovlp)[,1]]),
                           as.data.frame(gr.hg38[as.matrix(TreatedLike.ovlp)[,2]])) 

gr.TreatedLike.overlap$AMP <- sapply(gr.TreatedLike.overlap$CN, function(x) ifelse(x=="AMP",1,0))
gr.TreatedLike.overlap$DEL <- sapply(gr.TreatedLike.overlap$CN, function(x) ifelse(x=="DEL",1,0))

mat.TreatedLike.cn.amp <- acast(gene~sample, value.var = "AMP",fun.aggregate = sum,
                           data=gr.TreatedLike.overlap)
mat.TreatedLike.cn.amp[mat.TreatedLike.cn.amp > 1] <- 1
mat.TreatedLike.cn.del <- acast(gene~sample, value.var = "DEL",fun.aggregate = sum,
                           data=gr.TreatedLike.overlap)
mat.TreatedLike.cn.del[mat.TreatedLike.cn.del > 1] <- 1

# - Generate SNV and Indel matrix

# Load OAC SNV data 
oac_SNVs.nonsynonymous <- read.delim("SNVs.nonsynonymous.OAC.txt")
oac_SNVs.nonsynonymous$SampleID <- sub("_vs_.*", "", oac_SNVs.nonsynonymous$Sample)

# Load OAC Indel data
oac_indels <- read.delim("indels.OAC.txt")
oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample)  

# Extract only TreatedLike.Sig17 Samples
TreatedLike_snv <- oac_SNVs.nonsynonymous[which(oac_SNVs.nonsynonymous$SampleID %in% TreatedLikeSamples),]
TreatedLike_indel <- oac_indels[which(oac_indels$SampleID %in% TreatedLikeSamples),]

snvs.TreatedLike.selected <- unique(TreatedLike_snv[which(TreatedLike_snv$Gene %in% drivers),])
indels.TreatedLike.selected <- unique(TreatedLike_indel[which(TreatedLike_indel$Gene %in% drivers),])

mat.TreatedLike.snvs <- array(0,c(length(drivers),length(unique(c(snvs.TreatedLike.selected$SampleID,indels.TreatedLike.selected$SampleID,colnames(mat.TreatedLike.cn.amp))))))
rownames(mat.TreatedLike.snvs) <- drivers
colnames(mat.TreatedLike.snvs) <- unique(c(snvs.TreatedLike.selected$SampleID,indels.TreatedLike.selected$SampleID,colnames(mat.TreatedLike.cn.amp)))
for (i in 1:nrow(snvs.TreatedLike.selected)) {
  print(i)
  mat.TreatedLike.snvs[snvs.TreatedLike.selected[i,]$Gene,snvs.TreatedLike.selected[i,]$SampleID] <- 1
}

mat.TreatedLike.indels <- array(0,c(length(drivers),length(unique(c(snvs.TreatedLike.selected$SampleID,indels.TreatedLike.selected$SampleID,colnames(mat.TreatedLike.cn.amp))))))
rownames(mat.TreatedLike.indels) <- drivers
colnames(mat.TreatedLike.indels) <- unique(c(snvs.TreatedLike.selected$SampleID,indels.TreatedLike.selected$SampleID,colnames(mat.TreatedLike.cn.amp)))
for (i in 1:nrow(indels.TreatedLike.selected)) {
  print(i)
  mat.TreatedLike.indels[indels.TreatedLike.selected[i,]$Gene,indels.TreatedLike.selected[i,]$SampleID] <- 1
}


# - Generate Mat_list input for oncoprint heatmap

mat_list.TreatedLike <- list()
mat_list.TreatedLike$snv <- mat.TreatedLike.snvs[drivers,TreatedLikeSamples]
mat_list.TreatedLike$indel <- mat.TreatedLike.indels[drivers,TreatedLikeSamples]
mat.TreatedLike.cn.amp[drivers, TreatedLikeSamples][is.na(mat.TreatedLike.cn.amp[drivers, TreatedLikeSamples])] <- 0
mat_list.TreatedLike$amp <- mat.TreatedLike.cn.amp[drivers,TreatedLikeSamples]
mat.TreatedLike.cn.del[drivers, TreatedLikeSamples][is.na(mat.TreatedLike.cn.del[drivers, TreatedLikeSamples])] <- 0
mat_list.TreatedLike$del <- mat.TreatedLike.cn.del[drivers,TreatedLikeSamples]
names(mat_list.TreatedLike) <- c("SNV","Indel","AMP","DEL")

save(mat_list.TreatedLike, file = "mat_list.TreatedLike.Rdata")

# - Generate Oncoprint Heatmap

# Load Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0
mutsig.TreatedLike <- primaries_sigs_complete[TreatedLikeSamples,]

altered_nums<- Reduce( "+" , lapply(mat_list.TreatedLike, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.TreatedLike<-  lapply(mat_list.TreatedLike, function(x) x[slice_indx,])


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

pdf("oncoprint.drivers.TreatedLike.pdf",h=6,w=11)
oncoPrint(mat_list.TreatedLike,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS TreatedLike.Sig17+ Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.TreatedLike[,"SBS17a"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"))),  # Adjusted label size
            SBS17b = anno_barplot(mutsig.TreatedLike[,"SBS17b"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS2 = anno_barplot(mutsig.TreatedLike[,"SBS2"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS8 = anno_barplot(mutsig.TreatedLike[,"SBS8"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS18 = anno_barplot(mutsig.TreatedLike[,"SBS18"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS35 = anno_barplot(mutsig.TreatedLike[,"SBS35"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS40 = anno_barplot(mutsig.TreatedLike[,"SBS40"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS41 = anno_barplot(mutsig.TreatedLike[,"SBS41"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS44 = anno_barplot(mutsig.TreatedLike[,"SBS44"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          ))
dev.off()

