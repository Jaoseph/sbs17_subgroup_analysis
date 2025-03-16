# ----- Oncoprint for OCCAMS cohort EAC drivers

setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic Analysis/genomic_changes/Oncoprint")
library(GenomicRanges)
library(ComplexHeatmap)
library(reshape2)

drivers <- read.delim("OACdrivers_Frankell.txt", header=FALSE)$V1

# ----- NaiveLike.Sig17+ Samples

# Load sig17 status annotations, merge with mut sig contributions
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
NaiveLikeSamples <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"]

# - Generate CNV matrix

# Load OAC CNV data
oac_cna <- read.delim("CNAs.OAC.txt")
oac_cna$SampleID <- sub("_vs_.*", "", oac_cna$sample)

# Extract only NaiveLike samples
NaiveLike_cna <- oac_cna[which(oac_cna$SampleID %in% NaiveLikeSamples),]

# Load OAC driver gene locations
hg38 <- read.delim("hg38_geneLocations.txt", header = TRUE, sep = "\t")

gr.hg38 <- GRanges(
  seqnames = hg38$Chromosome.scaffold.name,
  ranges=IRanges(start=hg38$Gene.start..bp., end=hg38$Gene.end..bp.),
  gene=hg38$Gene.name)

gr.Naivelike.cns <- GRanges(
  seqnames = NaiveLike_cna$chromosome,
  ranges=IRanges(start=NaiveLike_cna$start, end=NaiveLike_cna$end),
  CN=NaiveLike_cna$CNchange,
  sample=NaiveLike_cna$SampleID)

# Find Gene corresponding to the CNA

NaiveLike.ovlp <- data.frame(findOverlaps(gr.Naivelike.cns,gr.hg38, type="any"))

gr.NaiveLike.overlap <- cbind(as.data.frame(gr.Naivelike.cns[as.matrix(NaiveLike.ovlp)[,1]]),
                           as.data.frame(gr.hg38[as.matrix(NaiveLike.ovlp)[,2]])) 

gr.NaiveLike.overlap$AMP <- sapply(gr.NaiveLike.overlap$CN, function(x) ifelse(x=="AMP",1,0))
gr.NaiveLike.overlap$DEL <- sapply(gr.NaiveLike.overlap$CN, function(x) ifelse(x=="DEL",1,0))

mat.NaiveLike.cn.amp <- acast(gene~sample, value.var = "AMP",fun.aggregate = sum,
                           data=gr.NaiveLike.overlap)
mat.NaiveLike.cn.amp[mat.NaiveLike.cn.amp > 1] <- 1
mat.NaiveLike.cn.del <- acast(gene~sample, value.var = "DEL",fun.aggregate = sum,
                           data=gr.NaiveLike.overlap)
mat.NaiveLike.cn.del[mat.NaiveLike.cn.del > 1] <- 1

# - Generate SNV and Indel matrix

# Load OAC SNV data 
oac_SNVs.nonsynonymous <- read.delim("SNVs.nonsynonymous.OAC.txt")
oac_SNVs.nonsynonymous$SampleID <- sub("_vs_.*", "", oac_SNVs.nonsynonymous$Sample)

# Load OAC Indel data
oac_indels <- read.delim("indels.OAC.txt")
oac_indels$SampleID <- sub("_vs_.*", "", oac_indels$Sample)  

# Extract only NaiveLike.Sig17 Samples
NaiveLike_snv <- oac_SNVs.nonsynonymous[which(oac_SNVs.nonsynonymous$SampleID %in% NaiveLikeSamples),]
NaiveLike_indel <- oac_indels[which(oac_indels$SampleID %in% NaiveLikeSamples),]

snvs.NaiveLike.selected <- unique(NaiveLike_snv[which(NaiveLike_snv$Gene %in% drivers),])
indels.NaiveLike.selected <- unique(NaiveLike_indel[which(NaiveLike_indel$Gene %in% drivers),])

mat.NaiveLike.snvs <- array(0,c(length(drivers),length(unique(c(snvs.NaiveLike.selected$SampleID,indels.NaiveLike.selected$SampleID,colnames(mat.NaiveLike.cn.amp))))))
rownames(mat.NaiveLike.snvs) <- drivers
colnames(mat.NaiveLike.snvs) <- unique(c(snvs.NaiveLike.selected$SampleID,indels.NaiveLike.selected$SampleID,colnames(mat.NaiveLike.cn.amp)))
for (i in 1:nrow(snvs.NaiveLike.selected)) {
  print(i)
  mat.NaiveLike.snvs[snvs.NaiveLike.selected[i,]$Gene,snvs.NaiveLike.selected[i,]$SampleID] <- 1
}

mat.NaiveLike.indels <- array(0,c(length(drivers),length(unique(c(snvs.NaiveLike.selected$SampleID,indels.NaiveLike.selected$SampleID,colnames(mat.NaiveLike.cn.amp))))))
rownames(mat.NaiveLike.indels) <- drivers
colnames(mat.NaiveLike.indels) <- unique(c(snvs.NaiveLike.selected$SampleID,indels.NaiveLike.selected$SampleID,colnames(mat.NaiveLike.cn.amp)))
for (i in 1:nrow(indels.NaiveLike.selected)) {
  print(i)
  mat.NaiveLike.indels[indels.NaiveLike.selected[i,]$Gene,indels.NaiveLike.selected[i,]$SampleID] <- 1
}


# - Generate Mat_list input for oncoprint heatmap

# Remove sample missing DEL and AMP data
NaiveLikeSamples <- NaiveLikeSamples[NaiveLikeSamples %in% colnames(mat.NaiveLike.cn.amp)]

mat_list.NaiveLike <- list()
mat_list.NaiveLike$snv <- mat.NaiveLike.snvs[drivers,NaiveLikeSamples]
mat_list.NaiveLike$indel <- mat.NaiveLike.indels[drivers,NaiveLikeSamples]
mat.NaiveLike.cn.amp[drivers, NaiveLikeSamples][is.na(mat.NaiveLike.cn.amp[drivers, NaiveLikeSamples])] <- 0
mat_list.NaiveLike$amp <- mat.NaiveLike.cn.amp[drivers,NaiveLikeSamples]
mat.NaiveLike.cn.del[drivers, NaiveLikeSamples][is.na(mat.NaiveLike.cn.del[drivers, NaiveLikeSamples])] <- 0
mat_list.NaiveLike$del <- mat.NaiveLike.cn.del[drivers,NaiveLikeSamples]
names(mat_list.NaiveLike) <- c("SNV","Indel","AMP","DEL")

save(mat_list.NaiveLike, file = "mat_list.NaiveLike.Rdata")

# - Generate Oncoprint Heatmap

# Load Mutational Signature Contributions
load("/Users/jao/Desktop/MSc_Project/OCCAMS_deconstructSigs/Primaries_deconstructSigs.Rdata")
primaries_sigs_complete[primaries_sigs_complete < 0.05] <- 0
mutsig.NaiveLike <- primaries_sigs_complete[NaiveLikeSamples,]

altered_nums<- Reduce( "+" , lapply(mat_list.NaiveLike, rowSums))

NUMBER_GENES<- 20

slice_indx<- order(altered_nums, decreasing = TRUE)[1:NUMBER_GENES]

mat_list.NaiveLike<-  lapply(mat_list.NaiveLike, function(x) x[slice_indx,])


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

pdf("oncoprint.drivers.NaiveLike.pdf",h=6, w=9)
oncoPrint(mat_list.NaiveLike,
          alter_fun = alter_fun_list, col = col,
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = "OCCAMS NaiveLike.Sig17+ Samples", column_title_side = "top",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          row_names_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = gpar(fontsize = 8, fontface = "bold"),
          top_annotation = HeatmapAnnotation(
            SBS17a = anno_barplot(mutsig.NaiveLike[,"SBS17a"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 6, fontface = "bold"))),  # Adjusted label size
            SBS17b = anno_barplot(mutsig.NaiveLike[,"SBS17b"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                  axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS2 = anno_barplot(mutsig.NaiveLike[,"SBS2"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS8 = anno_barplot(mutsig.NaiveLike[,"SBS8"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS18 = anno_barplot(mutsig.NaiveLike[,"SBS18"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS35 = anno_barplot(mutsig.NaiveLike[,"SBS35"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS40 = anno_barplot(mutsig.NaiveLike[,"SBS40"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS41 = anno_barplot(mutsig.NaiveLike[,"SBS41"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            SBS44 = anno_barplot(mutsig.NaiveLike[,"SBS44"], ylim = c(0, 0.7), height = unit(.7, "cm"), axis = TRUE, 
                                 axis_param = list(side = "left", gp = gpar(fontsize = 0))),
            annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
          ))
dev.off()

