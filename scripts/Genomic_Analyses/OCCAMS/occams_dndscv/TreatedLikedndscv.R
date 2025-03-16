# ------ Cohort-level dndscv - TreatedLike Samples
#setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")
library("dndscv")

load("TreatedLikegroup.Rdata")
dndsout.TreatedLike = dndscv(TreatedLike, outmats = T)
save(dndsout.TreatedLike, file = 'dndsout.TreatedLike.Rdata')

# dndscv outputs: Table of significant genes
sel_cv.TreatedLike = dndsout.TreatedLike$sel_cv

signif_genes.treated = sel_cv.TreatedLike[sel_cv.TreatedLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                       "wmis_cv","wnon_cv","wspl_cv")]

signif_genes.treated = sel_cv.TreatedLike[sel_cv.TreatedLike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv",
                                                                       "wmis_cv","wnon_cv","wspl_cv")]

