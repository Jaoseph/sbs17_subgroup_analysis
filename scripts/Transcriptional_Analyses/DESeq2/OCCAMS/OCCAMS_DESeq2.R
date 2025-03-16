##### DESeq2 differential expression analysis
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis")
library(ggrepel)

# Set Plot Theme
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# ----- Pre-processing for DESeq2 Analysis

# Load Raw expression counts

load("mat.expr.OCCAMS.RawCounts.RData")

# Load annotation status data

load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
all_samples <- occams_results_df$Sample
BELike <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"]
Naivelike <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"]
Treatedlike <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"]
Negative <- occams_results_df$Sample[occams_results_df$Phenotype_Assigned == "Sig17-"]

# Load sampleID annotations

load("rnaseq.plusonset.feb2024.RData")

# Find matching sequence_dna and sampleID

library(dplyr)

rnaseq.plusonset <- rnaseq.plusonset %>%
  filter(sequence_dna %in% all_samples) %>%
  select(SampleID, sequence_dna) #Total Samples with RNA-seq data: 272 ...

# Generate annotation object storing the matching sequence_dna and sampleID

annotations <- data.frame(
  SampleID = unique(rnaseq.plusonset$SampleID),
  sequence_dna = unique(rnaseq.plusonset$sequence_dna))

# Label the samples based on their Sig17 status

annotations <- annotations %>%
  mutate(Status = case_when(
    sequence_dna %in% BELike ~ "BarrettsLike",
    sequence_dna %in% Naivelike ~ "NaiveLike",
    sequence_dna %in% Treatedlike ~ "TreatedLike",
    sequence_dna %in% Negative ~ "Negative"))

table(annotations$Status)

# Extract samples with Status annotations 
mat.expr.OCCAMS <- as.data.frame(t(mat.expr.OCCAMS[rownames(mat.expr.OCCAMS) %in% annotations$SampleID,]))

# Only include protein-coding genes from GENCODE 

library(rtracklayer)

gencode <- import("gencode.v46.annotation.gtf", format = "gtf")
annotations_coding <- as.data.frame(mcols(gencode)) %>%
  filter(gene_type == "protein_coding")
coding_genes <- unique(annotations_coding$gene_name)
mat.expr.OCCAMS <- mat.expr.OCCAMS[rownames(mat.expr.OCCAMS) %in% coding_genes,] #counts input


# Check if the matrix and annotations are in the same order, rename Samples

all(colnames(mat.expr.OCCAMS) == annotations$SampleID)
colnames(mat.expr.OCCAMS) <- annotations$sequence_dna

# Remove SampleIDs from annotation and factorize status

sampleID <- annotations$sequence_dna
annotations <- factor(annotations[,"Status"], levels = c("NaiveLike", "BarrettsLike", "TreatedLike", "Negative")) #condition input
names(annotations) <- sampleID

# ----- DESeq2 Differential Expression Analysis

library(DESeq2)

# Re-format input into matrices and dataframe respectively

mat.expr.OCCAMS <- as.matrix(round(mat.expr.OCCAMS))
annotations <- as.data.frame(annotations)


# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = mat.expr.OCCAMS,
  colData = annotations,
  design = ~ annotations)

#keep <- rowSums(counts(dds)) >= 10
load("EdgeR_filterByExpr_genes.Rdata") #Ensures consistent genes filtered between DESeq2 and EdgeR
dds <- dds[keep,] #16223 genes 

# Run DESeq2 function
dds <- DESeq(dds)


# ----- Barretts vs Naive(Reference) Extract Results  ----- #

resLFC_BE <- results(dds, contrast = c("annotations", "BarrettsLike", "NaiveLike")) # logfc > 0 overexpresssed in BElike, logfc < 0 underexpressed in BElike
summary(resLFC_BE)

BE_dea <- as.data.frame(resLFC_BE)
BE_dea$Gene <- rownames(BE_dea); BE_dea <- BE_dea[,c(7,6,1,2,3,4,5)]; rownames(BE_dea) <- NULL

# Set FDR to 0.1 and l2fc to 0.5 and extract significant genes
fdrcutoff = 0.1
l2fccutoff = 0.5

write.table(as.data.frame(BE_dea[BE_dea$padj < fdrcutoff & abs(BE_dea$log2FoldChange) > l2fccutoff, ]),
            file="DESeq2/OCCAMS/BEvsNaive_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

## Volcano Plot
# Set expression color
BE_dea$diffexpressed <- "NO"
BE_dea$diffexpressed[BE_dea$log2FoldChange > 0.5 & BE_dea$padj < 0.1] <- "UP"
BE_dea$diffexpressed[BE_dea$log2FoldChange < -0.5 & BE_dea$padj < 0.1] <- "DOWN"

# Set gene labels
top25degs <- BE_dea[BE_dea$padj < 0.1 & abs(BE_dea$log2FoldChange) > 0.5, ]
top25degs <- head(top25degs[order(top25degs$padj), "Gene"], 25)
BE_dea$delabel <- ifelse(BE_dea$Gene %in% top25degs, BE_dea$Gene, NA)

# Generate Volcano Plot
pdf("DESeq2/OCCAMS/BEvsNaive_Volcano.pdf", width = 8, height = 9)
ggplot(data = BE_dea, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "grey", linetype = 'dashed') +
  geom_hline(yintercept = c(1), col = "grey", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  labs(color = 'BELike vs NaiveLike', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adjusted")) +
  geom_text_repel(max.overlaps = Inf)
dev.off()

# ----- Treated vs Naive(Reference) Extract Results  ----- #

resLFC_Treated <- results(dds, contrast = c("annotations", "TreatedLike", "NaiveLike")) # logfc > 0 overexpresssed in TreatedLike, logfc < 0 underexpressed in TreatedLike
summary(resLFC_Treated)

Treated_dea <- as.data.frame(resLFC_Treated)
Treated_dea$Gene <- rownames(Treated_dea); Treated_dea <- Treated_dea[,c(7,6,1,2,3,4,5)]; rownames(Treated_dea) <- NULL

Treated_dea <- na.omit(Treated_dea)

# Set FDR to 0.1 and l2fc to 0.5 and extract significant genes
fdrcutoff = 0.1
l2fccutoff = 0.5

write.table(as.data.frame(Treated_dea[Treated_dea$padj < fdrcutoff & abs(Treated_dea$log2FoldChange) > l2fccutoff, ]),
            file="DESeq2/OCCAMS/TreatedvsNaive_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

## Volcano Plot
# Set expression color
Treated_dea$diffexpressed <- "NO"
Treated_dea$diffexpressed[Treated_dea$log2FoldChange > 0.5 & Treated_dea$padj < 0.1] <- "UP"
Treated_dea$diffexpressed[Treated_dea$log2FoldChange < -0.5 & Treated_dea$padj < 0.1] <- "DOWN"

# Set gene labels
top25degs <- Treated_dea[Treated_dea$padj < 0.1 & abs(Treated_dea$log2FoldChange) > 0.5, ]
top25degs <- head(top25degs[order(top25degs$padj), "Gene"], 25)
Treated_dea$delabel <- ifelse(Treated_dea$Gene %in% top25degs, Treated_dea$Gene, NA)

# Generate Volcano Plot
pdf("DESeq2/OCCAMS/TreatedvsNaive_Volcano.pdf", width = 8, height = 9)
ggplot(data = Treated_dea, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "grey", linetype = 'dashed') +
  geom_hline(yintercept = c(1), col = "grey", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  labs(color = 'TreatedLike vs NaiveLike', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adjusted")) +
  geom_text_repel(max.overlaps = Inf)
dev.off()

# ----- Naive vs Sig17- Extract Results ----- #

resLFC_Negative <- results(dds, contrast = c("annotations", "NaiveLike", "Negative")) # logfc > 0 overexpresssed in Naivelike, logfc < 0 underexpressed in Naivelike
summary(resLFC_Negative)

Negative_dea <- as.data.frame(resLFC_Negative)
Negative_dea$Gene <- rownames(Negative_dea); Negative_dea <- Negative_dea[,c(7,6,1,2,3,4,5)]; rownames(Negative_dea) <- NULL

# Set FDR to 0.1 and l2fc to 0.5 and extract significant genes
fdrcutoff = 0.1
l2fccutoff = 0.5

write.table(as.data.frame(Negative_dea[Negative_dea$padj < fdrcutoff & abs(Negative_dea$log2FoldChange) > l2fccutoff, ]),
            file="DESeq2/OCCAMS/NaivevsNegative_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

write.table(Negative_dea, file="DESeq2/OCCAMS/NaivevsNegative_ALLGENES_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)


## Volcano Plot
# Set expression color
Negative_dea$diffexpressed <- "NO"
Negative_dea$diffexpressed[Negative_dea$log2FoldChange > 0.5 & Negative_dea$padj < 0.1] <- "UP"
Negative_dea$diffexpressed[Negative_dea$log2FoldChange < -0.5 & Negative_dea$padj < 0.1] <- "DOWN"

# Set gene labels
top25degs <- Negative_dea[Negative_dea$padj < 0.1 & abs(Negative_dea$log2FoldChange) > 0.5, ]
top25degs <- head(top25degs[order(top25degs$padj), "Gene"], 25)
Negative_dea$delabel <- ifelse(Negative_dea$Gene %in% top25degs, Negative_dea$Gene, NA)

# Generate Volcano Plot
pdf("DESeq2/OCCAMS/NaivevsNegative_Volcano.pdf", width = 8, height = 9)
ggplot(data = Negative_dea, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "grey", linetype = 'dashed') +
  geom_hline(yintercept = c(1), col = "grey", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  labs(color = 'NaiveLike vs Sig17-', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adjusted")) +
  geom_text_repel(max.overlaps = Inf)
dev.off()

# ----- Barretts vs Sig17- Extract Results ----- #

resLFC_BENegative <- results(dds, contrast = c("annotations", "BarrettsLike", "Negative")) # logfc > 0 overexpresssed in Naivelike, logfc < 0 underexpressed in Naivelike
summary(resLFC_BENegative)

BENegative_dea <- as.data.frame(resLFC_BENegative)
BENegative_dea$Gene <- rownames(BENegative_dea); BENegative_dea <- BENegative_dea[,c(7,6,1,2,3,4,5)]; rownames(BENegative_dea) <- NULL

# Set FDR to 0.1 and l2fc to 0.5 and extract significant genes
fdrcutoff = 0.1
l2fccutoff = 0.5

write.table(as.data.frame(BENegative_dea[BENegative_dea$padj < fdrcutoff & abs(BENegative_dea$log2FoldChange) > l2fccutoff, ]),
            file="DESeq2/OCCAMS/BarrettsvsNegative_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

write.table(BENegative_dea, file="DESeq2/OCCAMS/BarrettsvsNegative_ALLGENES_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

# Count up-regulated genes
num_upregulated <- nrow(BENegative_dea[BENegative_dea$padj < fdrcutoff & 
                                         BENegative_dea$log2FoldChange > l2fccutoff, ])

# Count down-regulated genes
num_downregulated <- nrow(BENegative_dea[BENegative_dea$padj < fdrcutoff & 
                                         BENegative_dea$log2FoldChange < -l2fccutoff, ])

# Print the results
cat("Number of significant up-regulated genes:", num_upregulated, "\n")
cat("Number of significant down-regulated genes:", num_downregulated, "\n")


## Volcano Plot
# Set expression color
BENegative_dea$diffexpressed <- "NO"
BENegative_dea$diffexpressed[BENegative_dea$log2FoldChange > 0.5 & BENegative_dea$padj < 0.1] <- "UP"
BENegative_dea$diffexpressed[BENegative_dea$log2FoldChange < -0.5 & BENegative_dea$padj < 0.1] <- "DOWN"

# Set gene labels
top25degs <- BENegative_dea[BENegative_dea$padj < 0.1 & abs(BENegative_dea$log2FoldChange) > 0.5, ]
top25degs <- head(top25degs[order(top25degs$padj), "Gene"], 25)
BENegative_dea$delabel <- ifelse(BENegative_dea$Gene %in% top25degs, BENegative_dea$Gene, NA)

# Generate Volcano Plot
pdf("DESeq2/OCCAMS/BarrettsvsNegative_Volcano.pdf", width = 8, height = 9)
ggplot(data = BENegative_dea, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "grey", linetype = 'dashed') +
  geom_hline(yintercept = c(1), col = "grey", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 8), xlim = c(-6, 6)) +
  scale_x_continuous(breaks = seq(-6, 6, 2)) +
  labs(color = 'BarrettsLike vs Sig17-', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adjusted")) +
  geom_text_repel(max.overlaps = Inf)
dev.off()


