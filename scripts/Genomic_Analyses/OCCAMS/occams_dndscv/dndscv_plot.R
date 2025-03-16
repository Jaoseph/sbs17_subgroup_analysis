setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")
library("dndscv")
library(readr)
library(ggpubr)
library(wesanderson)

# ------ Cohort-level dndscv - BELike Samples

# Run dndscv
#load("BELikegroup.Rdata")
#dndsout.BELike = dndscv(BELike, outmats = T)
#save(dndsout.BELike, file = 'dndsout.BELike.Rdata')

# Load dndscv output
load("dndsout.BELike.Rdata")

# dndscv outputs: Table of significant genes
sel_cv.BELike = dndsout.BELike$sel_cv

signif_genes.BELike = sel_cv.BELike[sel_cv.BELike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv","n_mis","n_non","n_spl","n_ind",
                                                                     "wmis_cv","wnon_cv","wspl_cv", "wind_cv")]

# calculate weighted dnds ratio
signif_genes.BELike <- signif_genes.BELike %>%
  mutate(weightedR = ((n_mis * wmis_cv) + (n_non * wnon_cv) + (n_spl * wspl_cv) + (n_ind * wind_cv))/(n_mis + n_non + n_spl + n_ind))
signif_genes.BELike <- subset(signif_genes.BELike, select = c(gene_name, qallsubs_cv, weightedR))
signif_genes.BELike$Phenotype <- "BELike"
  
# ------ Cohort-level dndscv - NaiveLike Samples

# Run dndscv
#load("NaiveLikegroup.Rdata")
#dndsout.NaiveLike = dndscv(NaiveLike, outmats = T)
#save(dndsout.NaiveLike, file = 'dndsout.NaiveLike.Rdata')

# Load dndscv output
load("dndsout.NaiveLike.Rdata")

# dndscv outputs: Table of significant genes
sel_cv.NaiveLike = dndsout.NaiveLike$sel_cv

signif_genes.naive = sel_cv.NaiveLike[sel_cv.NaiveLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv","n_mis","n_non","n_spl","n_ind",
                                                                          "wmis_cv","wnon_cv","wspl_cv", "wind_cv")]

# calculate weighted dnds ratio
signif_genes.naive <- signif_genes.naive %>%
  mutate(weightedR = ((n_mis * wmis_cv) + (n_non * wnon_cv) + (n_spl * wspl_cv) + (n_ind * wind_cv))/(n_mis + n_non + n_spl + n_ind))
signif_genes.naive <- subset(signif_genes.naive, select = c(gene_name, qallsubs_cv, weightedR))
signif_genes.naive$Phenotype <- "NaiveLike"

# ------ Cohort-level dndscv - TreatedLike Samples

# Run dndscv
#load("TreatedLikegroup.Rdata")
#dndsout.TreatedLike = dndscv(TreatedLike, outmats = T)
#save(dndsout.TreatedLike, file = 'dndsout.TreatedLike.Rdata')

# Load dndscv output
load("dndsout.TreatedLike.Rdata")

# dndscv outputs: Table of significant genes
sel_cv.TreatedLike = dndsout.TreatedLike$sel_cv

signif_genes.treated = sel_cv.TreatedLike[sel_cv.TreatedLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv","n_mis","n_non","n_spl","n_ind",
                                                                                 "wmis_cv","wnon_cv","wspl_cv", "wind_cv")]

# calculate weighted dnds ratio
signif_genes.treated <- signif_genes.treated %>%
  mutate(weightedR = ((n_mis * wmis_cv) + (n_non * wnon_cv) + (n_spl * wspl_cv) + (n_ind * wind_cv))/(n_mis + n_non + n_spl + n_ind))
signif_genes.treated <- subset(signif_genes.treated, select = c(gene_name, qallsubs_cv, weightedR))
signif_genes.treated$Phenotype <- "TreatedLike"

# ------ Cohort-level dndscv - Negative Samples

# Run dndscv
#load("negativegroup.Rdata")
#dndsout.negative = dndscv(negative, outmats = T)
#save(dndsout.negative, file = 'dndsout.negative.Rdata')

# Load dndscv output
load("dndsout.negative.Rdata")

# dndscv outputs: Table of significant genes
sel_cv.negative = dndsout.negative$sel_cv

signif_genes.negative = sel_cv.negative[sel_cv.negative$qallsubs_cv<0.1, c("gene_name","qallsubs_cv","n_mis","n_non","n_spl","n_ind",
                                                                          "wmis_cv","wnon_cv","wspl_cv", "wind_cv")]

# calculate weighted dnds ratio
signif_genes.negative <- signif_genes.negative %>%
  mutate(weightedR = ((n_mis * wmis_cv) + (n_non * wnon_cv) + (n_spl * wspl_cv) + (n_ind * wind_cv))/(n_mis + n_non + n_spl + n_ind))
signif_genes.negative <- subset(signif_genes.negative, select = c(gene_name, qallsubs_cv, weightedR))
signif_genes.negative$Phenotype <- "Negative"

# ------ Cohort-level dndscv Plot

# Combine all four data frames into one
combined_signif_genes <- rbind(signif_genes.BELike, 
                               signif_genes.naive, 
                               signif_genes.negative, 
                               signif_genes.treated)

combined_signif_genes$status <- ifelse(combined_signif_genes$Phenotype == "Negative", "Negative", "Positive")
combined_signif_genes$status <- factor(combined_signif_genes$status, levels = c("Positive", "Negative"))
combined_signif_genes$qallsubs_cv <- -log10(combined_signif_genes$qallsubs_cv) 
combined_signif_genes$qallsubs_cv[combined_signif_genes$qallsubs_cv == "Inf"] <- 15

my_colors <- c(
  "Negative" = wes_palette("AsteroidCity3")[2],
  "NaiveLike" = wes_palette("AsteroidCity3")[4],
  "BELike" = wes_palette("AsteroidCity3")[3],
  "TreatedLike" = wes_palette("AsteroidCity3")[1])

pdf("combined_dndscv.pdf", w = 7, h = 10)
ggscatter(combined_signif_genes, y = "qallsubs_cv", x = "weightedR",
          palette = my_colors,
          color = "Phenotype",
          label = combined_signif_genes$gene_name, repel = TRUE,
          facet.by = "status") + ylab("q-value")+
  xlab("Weighted dN/dS")
dev.off()



