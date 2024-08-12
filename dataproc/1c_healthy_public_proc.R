# Beta-value extraction for Otra et al. 2018 Healthy Breast tissue

rm(list=ls())
library(ENmix)
library(matrixStats)
library(minfi)

IDAT_DIR <- "*********** MASKED ***********"
OUTPUT_PATH <- "*********** MASKED ***********"
DET_P <- 0.000001; #ENmix default
SAMP_FRACTION <- 0.10
DPI <- 200

ChristensenLabMethArrayQCs <- function(rgSet, detP, sampThresh, dpi, outputPath, what2Return=c("ENmix","minfi"), usedFFPE=TRUE) {
  #'@description Phase 1 of 3: Run a set of existing MethylationEPIC QC pipelines
  if(getwd() != outputPath) setwd(outputPath)
  print(paste("Current working directory is:", getwd()))
  
  ## 1) Plot ENmix control plots:
  ENmix::plotCtrl(rgSet)
  
  ## 2) Run ENmix QC: 
  qc <- ENmix::QCinfo(rgSet, detPthre=detP, samplethre=sampThresh)
  
  ## 3) Generate minfi QC report:
  minfi::qcReport(rgSet=rgSet, pdf="minfi_QC_report.pdf")
  
  ## 4) minfi identification of outliers based on meth & unmeth intensities
  Mset <- minfi::preprocessRaw(rgSet)
  Mset <- minfi::minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE)
  png("minfi_OutliersByIntensity.png", height=8.27, width=11.69, units="in", res=dpi)
  minfi::plotQC(Mset$qc)
  title("Poor-performing Outlier Identification")
  dev.off()
  
  ## 6) Plot FFPE restoration probe:
  if(usedFFPE) {
    png("minfi_FFPEcontrol.png", height=8.27, width=11.69, units="in", res=dpi)
    controlStripPlot(rgSet, controls=c("RESTORATION"))
    dev.off()
  }
  
  if(what2Return == "ENmix") { 
    return(qc)
  } else if(what2Return == "minfi") {
    return(Mset)
  }
}

ChristensenLabMinfiAdoption <- function(rgSet, detP, sampThresh) {
  #'@description Phase 2 of 3: Executes minfi preprocessing by the Funnorm approach & excludes
  ## Step 1. Executes Funnorm  & methylumi.noob background correction:
  genomRatSet <- preprocessFunnorm(rgSet)
  print(genomRatSet)
  
  ## Step 2. Remove probes with low-call rate:
  print(paste(c("Detection P-value:","Sample threshold:"), c(detP, sampThresh))) 
  
  pvals <- detectionP(rgSet)
  failedP <- (pvals > detP) 
  
  print("Number of CpGs w/ failed P-values:")
  print(sum(failedP))
  
  print("Proportion of CpGs w/ failed P-values:")
  print(mean(failedP))
  
  print(paste("Number of CpGs w/ failed P-values in >", sampThresh, "of samples"))
  print(sum(rowMeans(failedP) > sampThresh))
  
  print(paste("Proportion of CpGs w/ failed P-values in >", sampThresh, "of samples"))
  print(mean(rowMeans(failedP) > sampThresh))
  
  failedProbes <- rownames(failedP)[rowMeans(failedP) > sampThresh]
  genomRatSet <- genomRatSet[! rownames(genomRatSet) %in% failedProbes]
  
  ## Step 3. Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
  genomRatSet <- dropMethylationLoci(genomRatSet) #drop technical probes
  genomRatSet <- dropLociWithSnps(genomRatSet) #drop SNPs w/ default MAF=0
  
  return(genomRatSet)
}

extractMethMatrix_final <- function(genomRatSet, targets=NULL, sampNameInTarget="full_filename", idToAdopt="Accession") {
  #'@description Extract beta-value matrix from preprocessed minfi objects
  methMat <- minfi::getBeta(genomRatSet)
  if(! is.null(targets) & identical(colnames(methMat),targets[,sampNameInTarget]) ) {
    print("Sample names matched! Proceed to barcode to ID switching...")
    colnames(methMat) <- targets$Accession
  }
  return(methMat);
}

main <- function() {
  print("************************ Process begins ************************")
  
  #------------------------------------Step 1. Load IDAT files as a RGChannelSet------------------------------------
  setwd(IDAT_DIR)
  targets <- read.metharray.sheet(getwd())
  targets$full_filename <- gsub("_Grn.idat", "", targets$full_filename)
  rgSet <- read.metharray.exp(targets=targets, extended=TRUE, force=TRUE)
  print(rgSet@annotation)
  
  #------------------------------------Step 2. Execute custom sets of QCs (minfi & ENmix)------------------------------------
  qc <- ChristensenLabMethArrayQCs(rgSet, detP=DET_P, sampThresh=SAMP_FRACTION, dpi=DPI, outputPath=OUTPUT_PATH, what2Return="ENmix", usedFFPE=FALSE)
  
  print("Total number of bad samples:")
  print(length(qc$badsample))
  
  ## Set aside probes to exclude:
  failedProbes <- qc$badCpG
  
  #------------------------------------Step 3. Basic raw-data exploration------------------------------------
  if(getwd() != OUTPUT_PATH) setwd(OUTPUT_PATH)
  
  ## Raw methylation density curves: 
  png("minfi_raw_densityCurves.png", height=8.27, width=11.69, units="in", res=DPI)
  par(mar=c(5,4,4,2))
  densityPlot(rgSet, main="Raw beta-value densities")
  dev.off()
  
  #------------------------------------Step 4. Run the customized minfi pipeline------------------------------------
  if(getwd() != OUTPUT_PATH) setwd(OUTPUT_PATH)
  genomRatSet <- ChristensenLabMinfiAdoption(rgSet, DET_P, SAMP_FRACTION)
  
  #------------------------------------Step 5. Plot final distribution & save------------------------------------
  OltraHealthyBetas <- extractMethMatrix_final(genomRatSet, targets=targets)
  save(list=c("OltraHealthyBetas","targets"), file="Oltra2018HealthySubjects.RData", compress=TRUE)
  
  ## MDS plots:
  png("minfi_mdsPlotPanels_post_normalization_1K_sites.png", height=8.27, width=11.69, units="in", res=DPI)
  mdsPlot(OltraHealthyBetas)
  dev.off()
  
  png("minfi_mdsPlotPanels_post_normalization_5K_sites.png", height=8.27, width=11.69, units="in", res=DPI)
  mdsPlot(OltraHealthyBetas, numPositions=5000)
  dev.off()
  
  ## Final methylation density curves:
  png("minfi_processed_density_curves.png", height=8.27, width=11.69, units="in", res=DPI)
  par(mar=c(5,4,4,2))
  densityPlot(OltraHealthyBetas, main="Preproocessed Methylation Densities")
  dev.off()
  
  print("************************ Process complete! ************************")
}

if(! interactive()) main()
