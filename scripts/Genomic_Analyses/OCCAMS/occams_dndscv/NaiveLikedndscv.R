# ------ Cohort-level dndscv - NaiveLike Samples
#setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")
library("dndscv")

load("NaiveLikegroup.Rdata")
dndsout.NaiveLike = dndscv(NaiveLike, outmats = T)
save(dndsout.NaiveLike, file = 'dndsout.NaiveLike.Rdata')

# dndscv outputs: Table of significant genes
sel_cv.NaiveLike = dndsout.NaiveLike$sel_cv

signif_genes.naive = sel_cv.NaiveLike[sel_cv.NaiveLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                       "wmis_cv","wnon_cv","wspl_cv")]

signif_genes.naive = sel_cv.NaiveLike[sel_cv.NaiveLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                       "wmis_cv","wnon_cv","wspl_cv")]

