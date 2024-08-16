# EWAS of TNBC Cohort 2 with Subject as Random Effect
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(limma)
source("../utils.R")

## Load & subset DNAm annotation & data:
SELE_UNIV <- read.table(PATH_UNIV)[,1] #25081 CpGs

ANNOT_EPIC <- loadEPICannotationFile()
ANNOT_EPIC <- subset(ANNOT_EPIC, Name %in% SELE_UNIV)

load_cohort2()
cohort2_betas <- subset(cohort2_betas, rownames(cohort2_betas) %in% SELE_UNIV)
cohort2_mvals <- minfi::logit2(cohort2_betas)
hist(cohort2_mvals) #verify M value distribution
stopifnot(identical(cohort2_covars$Tube.number, colnames(cohort2_mvals)))#checkpoint

## Design matrix:
myDesign <- model.matrix(~ MLPA.mean+AgeAtDx, data=cohort2_covars)
View(myDesign)

## Consensus correlation using **beta** values:
myBlock <- cohort2_covars$Patient.Number
dupCor <- duplicateCorrelation(cohort2_betas, myDesign, block=myBlock)
dupCor$consensus.correlation #0.687099
rm(cohort2_betas)

## EWAS w.r.t. **M** values:
fit <- lmFit(
  cohort2_mvals, 
  design = myDesign,
  correlation = dupCor$consensus.correlation, 
  block = myBlock
)

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

# write.csv(diff_loci, paste0(PATH_DMPS[["Cohort2"]]), row.names=FALSE, quote=FALSE)
