## Diffsig re-try. Sig17+ 5-FU vs Sig17+ Other treatment
setwd("/Users/jao/Desktop/MSc_Project/diffsig")

library(deconstructSigs)
library(diffsig, quietly = TRUE)
library(rstan, quietly = TRUE)
load("tutorial/fit.rda")
set.seed(241202)

# --- generate signature matrix object

# Load Signature Matrix
load("/Users/jao/Desktop/MSc_Project/Updated/sigs.allSamples.normalised.202101127_cosmicref13sig_S3S8.RData")
sigs.sbs <- colnames(sigs.allSamples) ; rm(sigs.allSamples)

# Load hg19 signatures
signatures.sbs96_hg19 <- read.table('COSMIC_v3.3.1_SBS_GRCh37.txt', h=T)
rownames(signatures.sbs96_hg19) <- signatures.sbs96_hg19$Type ; signatures.sbs96_hg19 <- signatures.sbs96_hg19[,-1]

C_k = signatures.sbs96_hg19[, names(signatures.sbs96_hg19) %in% sigs.sbs]

# --- Generate mutational counts (rows = Samples) and Risk Factor


# Load Clinical Data
clin <- read.csv("/Users/jao/Desktop/MSc_Project/Pan-cancer_Analysis/Analysis/Chemo_Analysis/Sigs_study_cohort_clinical demographics.csv")
clin$Sample <- sub("_vs_.*", "", clin$Sample)
# Extract only primaries
clin <- clin[clin$Category.Id == "PrimaryTumour",]
clin<- clin[,names(clin) %in% c("Sample", "Occams.ID", "Category.Id","Is.Chemo.Naive","TR.Chemotherapy.Treatment.Protocol")]
# Extract only treated samples
clin <- clin[clin$Is.Chemo.Naive == "false",]
# Remove Treated Samples Ambigious Treatments
clin[clin == ""] <- NA
clin <- clin[!is.na(clin$TR.Chemotherapy.Treatment.Protocol),]

FluoroTreat <- clin %>%
  filter(grepl("FU", clin$TR.Chemotherapy.Treatment.Protocol)) %>%
  pull(Sample) 

OtherTreat <- clin %>%
  filter(!grepl("FU", clin$TR.Chemotherapy.Treatment.Protocol)) %>%
  pull(Sample) 

Treated <- c(FluoroTreat, OtherTreat)

# Load OCCAMS cohort mutational counts
load("/Users/jao/Desktop/MSc_Project/Pan-cancer_Analysis/Analysis/oac_tally.Rdata")

# Read in stratified samples
#load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
#occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))

# Read in stratified samples
load("/Users/jao/Desktop/MSc_Project/Combined_Classifier/OCCAMS_CombinedresultsSummary.Rdata")
occams_results_df <- subset(occams_results_df, select = c(Sample,Phenotype_Assigned))
occams_results_df$Phenotype_Assigned <- ifelse(occams_results_df$Phenotype_Assigned == "Sig17-", "Sig17-", "Sig17+")

# Extract only the 721 primary samples and SBS contexts
oac_tally <- oac_tally[rownames(oac_tally) %in% occams_results_df$Sample,]
oac_tally <- oac_tally[,1:96]
oac_tally$Sample <- rownames(oac_tally)
oac_tally <- merge(occams_results_df, oac_tally, by = "Sample")

#Extract only Sig17+ Samples
#oac_tally <- oac_tally[oac_tally$Phenotype_Assigned == "TreatedLike.Sig17+",]
oac_tally <- oac_tally[oac_tally$Phenotype_Assigned == "Sig17+",]

oac_tally <- oac_tally[oac_tally$Sample %in% Treated,]
oac_tally$Treatment <- ifelse(oac_tally$Sample %in% FluoroTreat, 1, 0)
rownames(oac_tally) <- oac_tally$Sample
oac_tally <- oac_tally[,-c(1:2)]

# Generate Risk Factor X
X <- data.frame(
  intercept = rep(1, nrow(oac_tally)),
  treatment = oac_tally$Treatment
)
rownames(X) = rownames(oac_tally)
X <- as.data.frame(t(X))

# Generate Mutational Counts Y
Y = as.data.frame(t(oac_tally[,-97]))

# ----- Run DiffSig

identical(colnames(X), colnames(Y)) #check same sample order
identical(rownames(Y), rownames(C_k)) #check same mutational context order

fit <- diffsig_fit(X=X, Y=Y, C=C_k, beta_sd=0.5,
                   chains=4, cores=4, seed=241202)

# check Convergence
pdf("New_Retry/Positive_TreatmentComparisonTraceplot.pdf")
rstan::traceplot(fit, pars="beta")
dev.off()

pdf("New_Retry/Positive_TreatmentComparisonDiffSig.pdf")
diffsig_plot(fit, 
             pars=1:28, 
             ci_level = 80,
             signature_labels = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS8", "SBS17a", "SBS17b", "SBS18", "SBS28", "SBS30", "SBS35", "SBS40", "SBS41", "SBS44"), 
             riskfactor_labels = c("intercept","5-FU"),
             ncol=2, outlist = FALSE)
dev.off()



diffsig_plot <- function(fit, pars, signature_labels, riskfactor_labels,
                         est_color=NULL, colors=NULL,ncol=1,outlist=FALSE,ci_level=80) {
  
  ## Errors and warnings
  if(length(pars)!=length(signature_labels)*length(riskfactor_labels)) {
    stop("Number of parameters incorrect. length(pars) should equal to length(signature_labels)*length(riskfactor_labels).")
  }
  if(!is.null(est_color)) {
    if(length(est_color)!=1) {
      stop("Only 1 color can be specified for the estimation point est_color")
    }
    
  }
  if(!is.null(colors)) {
    if(length(unique(riskfactor_labels))!=length(colors)) {
      warning("Number of colors does not match the number of unique groups from riskfactor_labels")
    }
  }
  
  ## Setup data
  if(is.null(est_color)) {
    est_color <- "#FFC20A"
  }
  
  ci_low = (100-ci_level)/2/100
  ci_up = 1-ci_low
  
  statmat <- rstan::summary(fit, probs=c(ci_low,0.25,0.5,0.75,ci_up))$summary[,c("mean","10%", "25%", "50%", "75%", "90%")]
  statmat <- statmat[pars,]
  number_signatures <- length(signature_labels)
  number_riskfactors <- length(riskfactor_labels)
  statlist <- list()
  for (i in 1:number_riskfactors) {
    statlist[[i]] <- statmat[(number_signatures*i-(number_signatures-1)):(number_signatures*i),]
  }
  
  y <- as.numeric(seq(number_signatures, 1, by = -1))
  xlim.use <- c(min(statmat[, 2L]), max(statmat[, 6L]))
  xlim.use <- xlim.use + diff(xlim.use) * c(-0.05, 0.05)
  
  p.list <- list()
  for (i in 1:number_riskfactors) {
    xy.df <- data.frame(params = rownames(statlist[[i]]), y, statlist[[i]])
    xy.df$group <- rep(riskfactor_labels[[i]], times=number_signatures)
    colnames(xy.df) <- c("params", "y", "mean", "ll", "l", "m", "h", "hh","group")
    
    p.base <- ggplot2::ggplot(xy.df)
    p.name <- ggplot2::scale_y_continuous(breaks = y,
                                          labels = signature_labels,
                                          limits = c(0.8, y + 0.2))
    
    p.all <- p.base + ggplot2::xlim(xlim.use) + p.name + ggplot2::geom_vline(xintercept=0, linetype="dashed",color="darkgrey") +
      ggplot2::theme_bw()
    
    p.ci <- ggplot2::geom_segment(mapping = ggplot2::aes_string(x = "ll", xend = "hh", y = "y", yend = "y"))
    p.list[[i]] <- p.all + p.ci +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour="white", size=0.1),
                     panel.grid.minor.y = ggplot2::element_line(colour='grey', linetype='dashed', size=0.2))
  }
  
  if(!is.null(colors)) {
    if(length(unique(riskfactor_labels))!=length(colors)) {
      warning("Number of colors does not match the number of unique groups from riskfactor_labels")
    }
  }
  
  ## Colors not specified
  if(is.null(colors)) {
    color_by <- c("#929292","#2484F6","#00AD35","#E84A35","#CA001C","#C50077","#7C00FF")
    for(i in 1:number_riskfactors) {
      p.ci.2 <- ggplot2::geom_segment(ggplot2::aes_string(x = "l", xend = "h", y = "y", yend = "y"), color = color_by[i], size = 3)
      p.point <- ggplot2::geom_point(ggplot2::aes_string(x = "m", y = "y"), shape = 19, color=color_by[i], size = 3, show.legend = F)
      p.point2 <- ggplot2::geom_point(ggplot2::aes_string(x = "m", y = "y"), shape = 19, color=est_color, size = 2)
      
      p.list[[i]] <- p.list[[i]] +
        p.ci.2 +
        p.point +
        p.point2 +
        ggplot2::xlab(expression(beta)) +
        ggplot2::labs(caption=riskfactor_labels[i]) +
        ggplot2::theme(
          panel.grid.minor.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(angle = 0, hjust=1,margin = ggplot2::margin(t = -2, r = 0, b = 0, l = 0)),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(arrow = grid::arrow(length = ggplot2::unit(0.3, "cm"),
                                                                  ends = "both"), color="grey27"),
          axis.text.x = ggplot2::element_text(size="10"),
          panel.grid.major.x = ggplot2::element_line(),
          plot.margin = ggplot2::margin(t = 5,  # Top margin
                                        r = 10,  # Right margin
                                        b = 5,  # Bottom margin
                                        l = 5),
          panel.border = ggplot2::element_blank(),
          legend.margin = ggplot2::margin(t=-15),
          plot.caption= ggplot2::element_text(size=12, hjust=0.5, vjust=2, margin=ggplot2::margin(0,0,0,0))) +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
    }
  } else {
    if(length(colors)!=length(riskfactor_labels)) {
      stop("Length of colors should match the number of unique groups in riskfactor_labels")
    }
    
    for (i in 1:number_riskfactors) {
      p.ci.2 <- ggplot2::geom_segment(ggplot2::aes(x = l, xend = h, y = y, yend = y), color = colors[i], size = 3)
      p.point <- ggplot2::geom_point(ggplot2::aes(x = m, y = y), color=colors[i], shape = 19, size = 3, show.legend = F)
      p.point2 <- ggplot2::geom_point(ggplot2::aes_string(x = "m", y = "y"), color = est_color, shape = 19, size = 2)
      
      
      p.list[[i]] <- p.list[[i]] +
        p.ci.2 +
        p.point +
        p.point2 +
        ggplot2::xlab(expression(beta)) +
        ggplot2::labs(caption=riskfactor_labels[i]) +
        ggplot2::theme(
          panel.grid.minor.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(angle = 0, hjust=1,margin = ggplot2::margin(t = -2, r = 0, b = 0, l = 0)),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(arrow = grid::arrow(length = ggplot2::unit(0.3, "cm"),
                                                                  ends = "both"), color="grey27"),
          axis.text.x = ggplot2::element_text(size="10"),
          panel.grid.major.x = ggplot2::element_line(),
          plot.margin = ggplot2::margin(t = 5,  # Top margin
                                        r = 10,  # Right margin
                                        b = 5,  # Bottom margin
                                        l = 5),
          panel.border = ggplot2::element_blank(),
          legend.position = "none",
          plot.caption=ggplot2::element_text(size=12, hjust=0.5, vjust=2, margin=ggplot2::margin(0,0,0,0)))
    }
  }
  
  if (outlist==TRUE) {
    (p = p.list)
  } else {
    (p=cowplot::plot_grid(plotlist=p.list, ncol = ncol))
  }
  
  return(p)
}

