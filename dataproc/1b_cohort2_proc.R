# TNBC Cohort 2 methylation data processing
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(minfi)
library(ENmix)
library(matrixStats)

source("../utils.R") #procedures + imports constants/paths/gobalvar

## Script-specific constants:
SAMP_FRACTION <- 0.10 #custom
DET_P <- 0.000001 #ENmix default
DPI <- 200
QC_PATH <- "********** MASKED **********"

extractMethMatrix_final <- function(genomRatSet, targets) {
  methMat <- minfi::getBeta(genomRatSet)
  stopifnot(identical(colnames(methMat), targets$Complete_Barcode)) #checkpoint
  colnames(methMat) <- targets$Tube.number
  return(methMat)
}

Main <- function() {
  #------------------------------------ Load IDAT files as a RGChannelSet------------------------------------
  ## Load previously saved QCs, obtained by running the code once
  load(QC_PATH)
  rm(rgSet, targets) # will be reloaded

  ## Set aside samples to exclude:
  samplesToExclude <- qc$badsample
  samplesToExclude <- samplesToExclude[samplesToExclude != "202702320086_R05C01"] #decided to keep
  
  ## Set aside (initialize) probes to exclude:
  failedProbes <- qc$badCpG
  
  ## Reload IDATs EXCLUDING samples set aside:
  setwd(IDAT_DIR)
  targets <- read.metharray.sheet(getwd())
  targets <- subset(targets, ! Complete_Barcode %in% samplesToExclude)
  
  rgSet <- read.metharray.exp(targets=targets, extended=TRUE)
  
  #------------------------------------ minfi pipeline------------------------------------
  setwd(OUTPUT_PATH)
  ## Step 4.1. Executes Funnorm  & methylumi.noob background correction:
  genomRatSet <- preprocessFunnorm(rgSet)
  
  ## Step 4.2. Remove probes with low-call rate:
  print(paste(c("Detection P-value:","Sample threshold:"), c(DET_P, SAMP_FRACTION)))
  pvals <- detectionP(rgSet)
  failedP <- (pvals > DET_P)
  failedProbes <- union(failedProbes, rownames(failedP)[rowMeans(failedP) > SAMP_FRACTION]) #update
  
  ## Drop all failed CpGs detected by minfi: 
  print("Total number of failed probes:");
  print(length(failedProbes))
  genomRatSet <- genomRatSet[! rownames(genomRatSet) %in% failedProbes]
  
  ## Step 4.3. Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
  genomRatSet <- dropMethylationLoci(genomRatSet) #technical probes
  genomRatSet <- dropLociWithSnps(genomRatSet) #SNPs w/ default MAF=0
  print(genomRatSet); 
  
  #------------------------------------Save------------------------------------
  cohort2_betas <- extractMethMatrix_final(genomRatSet, targets=targets)
  save(list=c("cohort2_betas","targets"), file=OUT_PATH, compress=TRUE)
  
  ## Final methylation density plots:
  png("minfi_processed_density_curves.png", height=8.27, width=11.69, units="in", res=DPI)
  par(mar=c(5,4,4,2))
  densityPlot(cohort2_betas, main="Processed Methylation Densities", sampGroups=targets$MLPA)
  dev.off()
}

Main()
