# ------ Cohort-level dndscv - Negative Samples
#setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")
library("dndscv")

load("negativegroup.Rdata")
dndsout.negative = dndscv(negative, outmats = T)
save(dndsout.negative, file = 'dndsout.negative.Rdata')

# dndscv outputs: Table of significant genes
sel_cv.negative = dndsout.negative$sel_cv

signif_genes.negative = sel_cv.negative[sel_cv.negative$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                "wmis_cv","wnon_cv","wspl_cv")]

signif_genes.negative = sel_cv.negative[sel_cv.negative$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                           "wmis_cv","wnon_cv","wspl_cv")]

