# ------ Cohort-level dndscv - BELike Samples
setwd("/Users/jao/Desktop/MSc_Project/Combined_Classifier/Genomic_Analysis/dndscv")
library("dndscv")
library(readr)

# Run dndscv
load("BELikegroup.Rdata")
dndsout.BELike = dndscv(BELike, outmats = T)
save(dndsout.BELike, file = 'dndsout.BELike.Rdata')

# Run dndscv output
load("dndsout.BELike.Rdata")

# dndscv outputs: Table of significant genes
sel_cv.BELike = dndsout.BELike$sel_cv

signif_genes.BELike = sel_cv.BELike[sel_cv.BELike$qallsubs_cv<0.1, c("gene_name","qallsubs_cv","n_mis","n_non","n_spl","n_ind",
                                                                           "wmis_cv","wnon_cv","wspl_cv", "wind_cv")]

# calculate weighted dnds ratio

signif_genes.BELike <- signif_genes.BELike %>%
  mutate(weightedR = ((n_mis * wmis_cv) + (n_non * wnon_cv) + (n_spl * wspl_cv) + (n_ind * wind_cv))/(n_mis + n_non + n_spl + n_ind))




