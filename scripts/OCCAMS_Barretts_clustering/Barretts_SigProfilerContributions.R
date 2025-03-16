# Extract SBS and ID signatures in the ICGC cohort
setwd("/Users/jao/Desktop/MSc_Project/Barretts_Classifier")

library(tidyr)

# Load SBS data
sbs_pcawg <- read.csv('MutationalSignatureReferences/PCAWG_sigProfiler_SBS_signatures_in_samples.csv', row.names = 2)
sbs_nopcawg <- read.csv('MutationalSignatureReferences/nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv', row.names = 2)
sbs <- rbind(sbs_pcawg, sbs_nopcawg)


sbs_eac <- sbs[grepl(pattern = 'Eso', sbs$Cancer.Types), -c(1,2)]

sbs_eac_props <- data.frame(
  Sigs = colnames(sbs_eac),
  Prop_present = apply(sbs_eac, 2, function(x) sum(x > 0)/length(x))
)

# Load ID data
id <- read.csv('MutationalSignatureReferences/PCAWG_SigProfiler_ID_signatures_in_samples.csv', row.names = 2)
id_eac <- id[grepl(pattern = 'Eso', id$Cancer.Types), -c(1,2)]

id_eac_props <- data.frame(
  Sigs = colnames(id_eac),
  Prop_present = apply(id_eac, 2, function(x) sum(x > 0)/length(x))
)

# Collate and save
eac_sigs <- rbind(sbs_eac_props, id_eac_props)
write.table(eac_sigs, file = 'EAC_sigProfilerCont.txt',
            row.names = FALSE, col.names = TRUE)

