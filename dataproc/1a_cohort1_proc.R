# Pre-processing TNBC Cohort 1 Methylation EPIC 
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(minfi)
library(matrixStats)

source("../utils.R") #procedures + imports constants/paths/gobalvar

#------------------------ Load IDATs ------------------------
setwd(DIR_IDAT[["Cohort1"]])
targets <- read.metharray.sheet(getwd());
masterCovar <- read.csv(COVAR_DIR, stringsAsFactors=F)
colnames(masterCovar)[colnames(masterCovar) == "Complete.Barcode"] <- "Sample_Name"
masterCovar$Sample.ID <- NULL
targets <- merge(targets, masterCovar, by="Sample_Name")

## minfi preprocessing of tumors:
EPIC <- read.metharray.exp(targets=targets)
EPIC@annotation

## Extract sample ID and group label for plotting:
subjects <- pData(EPIC)$Sample_ID
groups <- pData(EPIC)$MLPA

## Density plots:
densityPlot(EPIC, main="Raw intensitie", sampGroups=groups, cex=0.75)
densityBeanPlot(EPIC, sampNames=subjects, sampGroups=groups, main="Raw intensities by sample")

## Outlier plot: Convert to a MethylSet
Mset <- preprocessRaw(EPIC)
Mset <- minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE)
plotQC(Mset$qc)
title("Poor-performing outlier identification")

#------------------------ Quality control & normalization ------------------------
## Step 1. Perform normalization (Funnorm) and background correction (methylumi.noob)
normalized850K <- preprocessFunnorm(EPIC)
dim(normalized850K)

## Step 2. Remove probes failed to meet detection P-value threshold of 0.05 in >20% samples:
pvals <- detectionP(EPIC)
failedP <- (pvals > 0.05) #Be careful with the inequality sign!

fraction <- 0.20
failedProbes <- rownames(failedP)[rowMeans(failedP) > fraction] #list of probes
sum(rowMeans(failedP) > fraction)

normalized850K <- normalized850K[! rownames(normalized850K) %in% failedProbes]
dim(normalized850K)

## Step 3. Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
normalized850K <- dropMethylationLoci(normalized850K) #drop technical SNP probes
normalized850K <- dropLociWithSnps(normalized850K) #drop polymorphic SNPs; MAF set to 0 by default
dim(normalized850K)

## Step 4. Extract beta values & switch out subject ID
EPIC.betas <- getBeta(normalized850K)
EPIC.betas <- EPIC.betas[ , match(targets$Sample_Name, colnames(EPIC.betas))]
colnames(EPIC.betas) <- as.character(targets$Sample_ID)
stopifnot(identical(as.character(targets$Sample_ID), colnames(EPIC.betas)))

## Visualize densities of the final data set:
densityPlot(EPIC.betas, sampGroups=groups, main="Normalized, low-call-rate/SNP/non-CpG sites removed");

## Export:
save(
  list = c("EPIC.betas", "targets"),
  file = BETAS_PATH,
  compress = TRUE
)
