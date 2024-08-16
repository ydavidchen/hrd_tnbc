# RAD51C array CpG Visualization
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(pheatmap)
source("../utils.R")
source("../plotThemes.R")

TSS_RAD51C <- 56770931 #hg19
BP <- 1500

## Load methylation data & clinical:
load_cohort1()
load_cohort2()
rownames(cohort1_covars) <- cohort1_covars$Sample_ID
rownames(cohort2_covars) <- cohort2_covars$Tube.number

## Load MethylationEPIC Annotation & Identify PolakCpG CpGs:
ANNOT_RAD51C <- loadEPICannotationFile()
ANNOT_RAD51C <- subset(ANNOT_RAD51C, grepl("RAD51C", UCSC_RefGene_Name))
ANNOT_RAD51C$PolakCpG <- ANNOT_RAD51C$pos >= TSS_RAD51C-BP & ANNOT_RAD51C$pos <= TSS_RAD51C+BP #upstream & downstream
sum(ANNOT_RAD51C$PolakCpG) #16

## Subset & Order by genomic position:
common_cpgs <- Reduce(intersect, list(ANNOT_RAD51C$Name, rownames(cohort1_betas), rownames(cohort2_betas)))
length(common_cpgs)

ANNOT_RAD51C <- ANNOT_RAD51C[common_cpgs, ]
ANNOT_RAD51C <- ANNOT_RAD51C[order(ANNOT_RAD51C$pos), ]

cohort1_betas <- cohort1_betas[ANNOT_RAD51C$Name, ]
cohort2_betas <- cohort2_betas[ANNOT_RAD51C$Name, ]

stopifnot( identical(rownames(cohort1_betas), ANNOT_RAD51C$Name) ) #checkpoint
stopifnot( identical(rownames(cohort2_betas), ANNOT_RAD51C$Name) ) #checkpoint

## Heatmap plotting objects:
CPG_ANNOT <- createCpGTrackingBars(ANNOT_RAD51C)
CPG_ANNOT$PolakCpG <- ifelse(rownames(CPG_ANNOT) %in% ANNOT_RAD51C$Name[ANNOT_RAD51C$PolakCpG], "Yes", "No")

ANNOT_COLORS <- list(Context = c(Island="black", Shore="dimgray", Shelf="gray60", OpenSea="gray90"))
ANNOT_COLORS[["HRD"]] <- ANNOT_COLORS[["PolakCpG"]] <- COLORS_BINARY

## Heat map:
hm_wrapper <- function(mat, samp_annot, title) {
  pheatmap(
    mat,
    cutree_cols = 2,
    cluster_rows = FALSE, #keeps manual order
    annotation_col = samp_annot,
    annotation_row = CPG_ANNOT[ , c("Context","PolakCpG")],
    annotation_colors = ANNOT_COLORS,
    clustering_distance_rows = CL_PARAMS[1],
    clustering_distance_cols = CL_PARAMS[1],
    clustering_method = CL_PARAMS[2], 
    color = HEAT_COLS,
    border_color = NA,
    fontsize = 12,
    fontsize_row = 7, 
    main = title
  )
}

CL_PARAMS <- c("euclidean", "ward.D2")
hm_wrapper(cohort1_betas, cohort1_covars[,"HRD",drop=F], "TNBC Cohort 1")
hm_wrapper(cohort2_betas, cohort2_covars[,"HRD",drop=F], "TNBC Cohort 2")
