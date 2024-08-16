# EWAS of TNBC Cohort 1 (all distinct subjects)
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(limma)
source("../utils.R")

## Load & subset data & annotations:
SELE_UNIV <- read.table(PATH_UNIV)[,1] #25081 CpGs

ANNOT_EPIC <- loadEPICannotationFile()
ANNOT_EPIC <- subset(ANNOT_EPIC, Name %in% SELE_UNIV)

load_cohort1()
cohort1_betas <- subset(cohort1_betas, rownames(cohort1_betas) %in% SELE_UNIV)
cohort1_mvals <- minfi::logit2(cohort1_betas)
hist(cohort1_mvals) #verify distribution
stopifnot(identical(cohort1_covars$Sample_ID, colnames(cohort1_mvals))) #checkpoint

## Design matrix:
myDesign <- model.matrix( ~ MLPA.mean+AgeAtDx, data=cohort1_covars)
View(myDesign)

## EWAS w.r.t. **M** values:
fit <- lmFit(cohort1_mvals, design=myDesign)
fit <- eBayes(fit)

## Extract & export results:
diff_loci <- topTable(
  fit,
  number = Inf,
  coef = "MLPA.mean", 
  genelist = ANNOT_EPIC, 
  adjust.method = "fdr", 
  sort.by = "p"
)

# write.csv(diff_loci, paste0(PATH_DMPS[["Cohort1"]]), row.names=FALSE, quote=FALSE)
