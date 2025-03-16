##### DESeq2 differential expression analysis
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Transcriptional_Analysis")
library(TCGAbiolinks)
library(ggrepel)
library(rtracklayer)

# Set Plot Theme
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)))

# Load TCGA Sig17 classification annotation
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/CombinedClassifer_TCGA_resultsSummary.Rdata")

# ----- Pre-processing for DESeq2 Analysis

# Load TCGA Raw expression counts 

#Query
query <- GDCquery(
  project = 'TCGA-ESCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = results.tcga_df$Patient
)
#GDCdownload(query)
expr.tcga <- GDCprepare(query = query)
expr.tcga <- expr.tcga[,expr.tcga$sample_type == 'Primary Tumor'] #only primary tumors
expr.tcga <- expr.tcga[!duplicated(rowData(expr.tcga)$gene_name) &
                         !is.na(rowData(expr.tcga)$gene_name), ] #remove duplicates

#Extract raw counts
expr.tumor <- assay(expr.tcga, 'unstranded')
rownames(expr.tumor) <- rowData(expr.tcga)$gene_name #Re-format names
colnames(expr.tumor) <- sapply(colnames(expr.tumor),
                               function(x) substr(x,1,12)) #Re-format sampleIDs

#Only include protein-coding genes from GENCODE 
gencode <- import("/Users/jao/Desktop/MSc_Project/Pan-cancer_Analysis/Analysis/gencode.v46.annotation.gtf", format = "gtf")
annotations_coding <- as.data.frame(mcols(gencode)) %>%
  filter(gene_type == "protein_coding")
coding_genes <- unique(annotations_coding$gene_name)
expr.tumor <- expr.tumor[rownames(expr.tumor) %in% coding_genes,] #19116 coding genes

# Generate sample annotations for DESeq2 input
annotations <- data.frame(
  Status = results.tcga_df$Phenotype_Assigned,
  row.names = results.tcga_df$Patient
)

annotations$Status <- ifelse(annotations$Status == "NaiveLike.Sig17+", "NaiveLike", "Negative")

annotations <- factor(annotations$Status, levels = c("Negative","NaiveLike")) #factorize
names(annotations) <- results.tcga_df$Patient 

# Check if the matrix and annotations are in the same order
all(rownames(annotations) == colnames(expr.tumor))

# ----- DESeq2 Differential Expression Analysis

library(DESeq2)

# Re-format input into matrices and dataframe respectively

expr.tumor <- as.matrix(round(expr.tumor))
annotations <- as.data.frame(annotations)

# Generate DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = expr.tumor,
  colData = annotations,
  design = ~ annotations
)
#keep <- rowSums(counts(dds)) >= 10
load('tcga_EdgeR_filterByExpr_genes.Rdata') #Ensures consistent genes filtered between DESeq2 and EdgeR
dds <- dds[keep,] #16642

# Run DESeq2 function
dds <- DESeq(dds)

# Extract Results
resultsNames(dds)
res <- results(dds, name = "annotations_NaiveLike_vs_Negative") # logfc > 0 overexpresssed in Naivelike, logfc < 0 underexpressed in Naivelike
summary(res)

res <- as.data.frame(res)
res$Gene <- rownames(res); res <- res[,c(7,6,1,2,3,4,5)]; rownames(res) <- NULL

# Set FDR to 0.1 and l2fc to 0.5 and extract significant genes
fdrcutoff = 0.1
l2fccutoff = 0.5

write.table(as.data.frame(res[res$padj < fdrcutoff & abs(res$log2FoldChange) > l2fccutoff, ]),
            file="DESeq2/TCGA/NaivevsNegative_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

write.table(res, file="DESeq2/TCGA/NaivevsNegative_ALLGENES_DESeq2_output.tsv", sep="\t", quote=F, col.names = T)

## Volcano Plot
# Set expression color
res$diffexpressed <- "NO"
res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.1] <- "UP"
res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.1] <- "DOWN"


# Set gene labels
top25degs <- res[res$padj < 0.1 & abs(res$log2FoldChange) > 0.5, ]
top25degs <- head(top25degs[order(top25degs$padj), "Gene"], 25)
res$delabel <- ifelse(res$Gene %in% top25degs, res$Gene, NA)

# Generate Volcano Plot
pdf("DESeq2/TCGA/NaivevsNegative_Volcano.pdf", width = 8, height = 9)
ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
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



