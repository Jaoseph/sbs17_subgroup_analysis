# Compare Thymidylate Synthetase (TYMS) expression between NaiveLike and TreatedLike?
library(dplyr)

##### Transcriptional Analysis pre-processing (TMM, CPM, Log transformation)

setwd("/Users/jao/Desktop/MSc_Project/CombinedClassifier_ExprAnalysis")

# Load Raw expression counts
load("mat.expr.OCCAMS.RawCounts.RData")

# Load Samples 
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
primaries <- occams_results_df$Sample
rownames(occams_results_df) <- occams_results_df$Sample; occams_results_df <- subset(occams_results_df, select = Phenotype_Assigned)

# Extract matching sequence_dna and sampleID
load("rnaseq.plusonset.feb2024.RData")
rnaseq.plusonset <- rnaseq.plusonset %>%
  filter(sequence_dna %in% primaries) %>%
  select(SampleID, sequence_dna) #272 samples

# Generate annotation object storing the matching sequence_dna and sampleID
annotations <- data.frame(
  SampleID = unique(rnaseq.plusonset$SampleID),
  sequence_dna = unique(rnaseq.plusonset$sequence_dna))

# Add Sig17 Status
annotations$Status <- ifelse(annotations$sequence_dna %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "NaiveLike.Sig17+"], "NaiveLike.Sig17+",
                                 ifelse(annotations$sequence_dna %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "BarrettsLike.Sig17+"], "BarrettsLike.Sig17+",
                                        ifelse(annotations$sequence_dna %in% rownames(occams_results_df)[occams_results_df$Phenotype_Assigned == "TreatedLike.Sig17+"], "TreatedLike.Sig17+", "Sig17-")))


# Exract raw expression counts for primaries
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
annotations <- factor(annotations[,"Status"]) #condition input
names(annotations) <- sampleID

# ----- Differential expression / Transcriptional analysis EdgeR Pre-processing 

suppressWarnings(library(edgeR, quietly = T))

# Generate an DGEList object 

y <- DGEList(counts = mat.expr.OCCAMS, group = annotations)

# Remove genes with consistently very low or zero counts

keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]

# Perform TMM normalization and transfer to CPM (Counts per Million)

y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y, log = TRUE, prior.count = 1)
count_norm<-as.data.frame(count_norm) #15299 coding genes
filterByExpr_genes <- rownames(count_norm)

# Save

save(count_norm, file = "Processed_ExpressionCounts.Rdata")
save(annotations, file = "Processed_ExpressionCounts_Annotations.Rdata")
save(keep, file = "EdgeR_filterByExpr_genes.Rdata") #maintain same number of genes for both ML and DEA

annotations <- data.frame(
  status = annotations,
  row.names = names(annotations)
)
annotations <- as.data.frame(t(annotations))

expr.matrix <- as.data.frame(t(rbind(annotations, count_norm)))
save(expr.matrix, file = "Processed_ExprMatrix.Rdata")
