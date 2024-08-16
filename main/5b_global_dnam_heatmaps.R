# Heatmaps for Global Methylation
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(pheatmap)
source("../utils.R")
source("../plotThemes.R")

SELE_UNIV <- read.table(PATH_UNIV)[,1]

CPG_ANNOT <- createCpGTrackingBars()

ANNOT_COLORS[["BRCA1"]] <- c(Hyperme="black", Mutated="lightgray") #from plotThemes

wrapper_hm <- function(mat, samp_annot, gaps_col=NULL) {
  #'@description Heatmap w/ optional custom ordering
  #'@param gaps_col Position to make gaps. If NULL, order ignored & hclust will be run
  ph <- pheatmap(
    mat,
    show_rownames = FALSE,
    cluster_cols = is.null(gaps_col),
    treeheight_row = 0, #hide
    border_color = NA,
    fontsize = 12,
    color = HEAT_COLS,
    annotation_colors = ANNOT_COLORS,
    annotation_row = CPG_ANNOT,
    annotation_col = samp_annot,
    gaps_col = gaps_col
  )
}

# ---------------------- Load & Re-order/Match Cohort Data ---------------------- 
load_cohort1()
load_cohort2()

rownames(cohort1_covars) <- cohort1_covars$Sample_ID
rownames(cohort2_covars) <- cohort2_covars$Tube.number

cohort1_covars <- cohort1_covars[order(cohort1_covars$RPMMSampleOrder), ]
cohort2_covars <- cohort2_covars[order(cohort2_covars$RPMMSampleOrder), ]

cohort1_betas <- cohort1_betas[SELE_UNIV, match(rownames(cohort1_covars), colnames(cohort1_betas))]
cohort2_betas <- cohort2_betas[SELE_UNIV, match(rownames(cohort2_covars), colnames(cohort2_betas))]

# ---------------------- Heat Maps w/ Unsupervised Clustering ----------------------
## Global methylation w/ RPMM sample order across all samples (Fig.2): 
wrapper_hm(cohort1_betas,  cohort1_covars[,c("HRD","RPMM")],  c(17,25))
wrapper_hm(cohort2_betas,  cohort2_covars[,c("HRD","RPMM")],  c(35,48))


## EWAS loci w/ hclust (Suppl. Fig):
EWAS_CONCORD <- read.csv(PATH_DMPS[["both"]])

wrapper_hm(cohort1_betas[EWAS_CONCORD$Name,], cohort1_covars[,c("HRD","AgeAtDx")])
wrapper_hm(cohort2_betas[EWAS_CONCORD$Name,], cohort2_covars[,c("HRD","AgeAtDx")])


## Global methylation for HRD BRCA1-altered w/ hclust (Reviewer Request):
cohort1_sele <- subset(cohort1_covars, HRD=="Yes" & (BRCA.Status=="BRCA1" | BRCA1_Hypermeth=="Yes"))
cohort2_sele <- subset(cohort2_covars, HRD=="Yes" & (BRCA.Status=="BRCA1" | isBRCA1Meth))

cohort1_sele$BRCA1 <- ifelse(cohort1_sele$BRCA1_Hypermeth=="Yes", "Hyperme", "Mutated")
cohort2_sele$BRCA1 <- ifelse(cohort2_sele$isBRCA1Meth, "Hyperme", "Mutated")

wrapper_hm(cohort1_betas[,rownames(cohort1_sele)], cohort1_sele[,"BRCA1",drop=F])
wrapper_hm(cohort2_betas[,rownames(cohort2_sele)], cohort2_sele[,"BRCA1",drop=F])
