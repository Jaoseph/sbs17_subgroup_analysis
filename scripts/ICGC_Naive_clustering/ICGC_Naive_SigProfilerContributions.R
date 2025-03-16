# Extract SBS and ID signatures in the ICGC cohort
setwd("/Users/jao/Desktop/MSc_Project/ChemoNaive_Only")

library(tidyr)

# Load SBS data
sbs_pcawg <- read.csv('MutationalSignatureReferences/PCAWG_sigProfiler_SBS_signatures_in_samples.csv', row.names = 2)
sbs_nopcawg <- read.csv('MutationalSignatureReferences/nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv', row.names = 2)
sbs <- rbind(sbs_pcawg, sbs_nopcawg)

icgc_wgs_cancer_types <- c("Biliary-AdenoCA", "Biliary-AdenoCa", "Breast-AdenoCA", "Breast-cancer",
                           "Breast-DCIS", "Breast-LobularCA", "Eso-AdenoCA",
                           "Head-SCC", "Eso-AdenoCa", "Kidney-RCC", "Liver-HCC",
                           "Lymph-BNHL", "Myeloid-AML", "AML", "Ovary-AdenoCA", "Ovary-AdenoCa",
                           "Panc-AdenoCA", "Panc-AdenoCa", "Prost-AdenoCA", "Prost-AdenoCa",
                           "Skin-Melanoma", "Stomach-AdenoCA", "Stomach-AdenoCa")



sbs_icgc <- sbs[sbs$Cancer.Types %in% icgc_wgs_cancer_types, -c(1,2)]

sbs_icgc_props <- data.frame(
  Sigs = colnames(sbs_icgc),
  Prop_present = apply(sbs_icgc, 2, function(x) sum(x > 0)/length(x))
)

# Load ID data
id <- read.csv('MutationalSignatureReferences/PCAWG_SigProfiler_ID_signatures_in_samples.csv', row.names = 2)
id_icgc <- id[id$Cancer.Types %in% icgc_wgs_cancer_types, -c(1,2)]

id_brca_props <- data.frame(
  Sigs = colnames(id_icgc),
  Prop_present = apply(id_icgc, 2, function(x) sum(x > 0)/length(x))
)

# Collate and save
icgc_sigs <- rbind(sbs_icgc_props, id_brca_props)
write.table(icgc_sigs, file = 'NaiveOnly_sigProfilerCont.txt',
            row.names = FALSE, col.names = TRUE)

