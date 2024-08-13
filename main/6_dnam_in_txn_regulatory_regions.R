# Methylation around Transcriptional Regulatory Regions

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(matrixStats)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
source("../utils.R")
source("../plotThemes.R")

RANGE_SIZE <- 1000 #bp
CHROMS <- paste0("chr", 1:22)
SMOOTH_FUNC <- function(x) stats::lowess(x, f=1/5) #lambda func, lowess not loess

## Paths to ENCODE reference profiles:
HMEC_H3K27AC <- paste0(DIR_ENCODE, "GSM733660_hg19_wgEncodeBroadHistoneHmecH3k27acStdPk.broadPeak.gz")
HMEC_H3K36ME3 <- paste0(DIR_ENCODE, "GSM733707_hg19_wgEncodeBroadHistoneHmecH3k36me3StdPk.broadPeak.gz")

## GRanges/Genomation wrappers:
restrict_on_seqnames <- function(gr, sn2subset=CHROMS) {
  #'@description Helper function to clean GRanges object
  gr <- gr[gr@seqnames %in% sn2subset, ]
  seqlevels(gr) <- seqlevelsInUse(gr)
  return(gr)
}

proc_filter_encode_chip <- function(path, adjMeth="bonferroni", qThresh=0.05) {
  #'@description Helper function to load ENCODE data
  #'@describeIn https://genome.ucsc.edu/FAQ/FAQformat.html
  broadPeak <- rtracklayer::import(path)
  broadPeak <- restrict_on_seqnames(broadPeak)
  
  broadPeak$pValue <- 10^(-broadPeak$pValue);
  broadPeak$qValue <- p.adjust(broadPeak$pValue, method=adjMeth)
  broadPeak <- broadPeak[broadPeak$qValue < qThresh, ]
  return(granges(broadPeak, use.mcols=FALSE))
}

assemble_summ_beta_granges <- function(mat, annot_genomation, agg_func=rowMedians) {
  #'@description Assembles GRanges object
  #'@param mat Matrix of beta-values with row.names=cgID, col.names=samples
  #'@param annot_genomation Illumina annotation with Name (CpG name), chr, start, end, columns
  #'@param agg_func 2D aggregation functions, e.g. rowMeans or rowMedians
  betas_summary <- data.frame(Name=row.names(mat), beta=agg_func(mat)) 
  betas_summary <- merge(annot_genomation, betas_summary, by="Name")
  gr_obj <- makeGRangesFromDataFrame(betas_summary, keep.extra.columns=TRUE)
  gr_obj$Name <- NULL
  return(gr_obj)
}

# ------------------------------ Build Genomation Reference Windows ------------------------------
## Illumina CpG annotation:
ANNOT_850K <- loadEPICannotationFile()
annot_genomation <- data.frame(
  Name = ANNOT_850K$Name,
  chrom = ANNOT_850K$chr,
  strand = ANNOT_850K$strand,
  start = as.integer(ANNOT_850K$pos)
)
annot_genomation$end <- annot_genomation$start

## Annotated hg19 TSS from UCSC, avail. in Bioconductor:
TRANSCRIPT_REGIONS <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
TRANSCRIPT_REGIONS <- restrict_on_seqnames(TRANSCRIPT_REGIONS)
TRANSCRIPT_REGIONS$tx_id <- TRANSCRIPT_REGIONS$tx_name <- NULL
end(TRANSCRIPT_REGIONS) <- start(TRANSCRIPT_REGIONS)
TRANSCRIPT_REGIONS <- flank(TRANSCRIPT_REGIONS, width=RANGE_SIZE, both=TRUE) #strand already taken into account
TRANSCRIPT_REGIONS
range(width(TRANSCRIPT_REGIONS))

## HMEC ChIP-seq annotations from ENCODE (see data dictionary):
PEAKS_HMEC_H3K27ac <- proc_filter_encode_chip(HMEC_H3K27AC)
PEAKS_HMEC_H3K27ac
summary(width(PEAKS_HMEC_H3K27ac))

PEAKS_HMEC_H3K36me3 <- proc_filter_encode_chip(HMEC_H3K36ME3)
PEAKS_HMEC_H3K36me3
summary(width(PEAKS_HMEC_H3K36me3))

# --------------------------------------- Execution ---------------------------------------
wrapper_genom <- function(gr1, gr0, txdb_annot, bin_num=100,
                          window_name="Distance to Genomic Feature (bp)",
                          gr_title=c("HRD","Non-HRD"), cols=c("red","blue"),
                          ranges=c(-RANGE_SIZE,RANGE_SIZE), meth_var_name="Summary beta-value", pltTitle="") {
  #'@description Sketches 2 overlapping curves of summarized beta-values
  #'@param g1,g0 GRanges object
  #'@param txdb_annot GenomicFeatures extracted from TxDb objects
  
  ## Plotting data object:
  gr_list <- list(gr1, gr0)
  names(gr_list) <- gr_title
  print(names(gr_list))
  
  sml <- ScoreMatrixList(
    targets = gr_list,
    windows = txdb_annot, 
    bin.num = bin_num,
    bin.op = "median",
    weight.col = "beta",
    is.noCovNA = TRUE,
    strand.aware = TRUE
  )
  
  ## Overalpping curves in 1 figure:
  plotMeta(
    sml,
    smoothfun = SMOOTH_FUNC,
    xcoords = ranges,
    line.col = cols,
    main = pltTitle,
    xlab = window_name, 
    ylab = meth_var_name, 
    lwd = 3, 
    bty = "L", 
    cex = 1.5, 
    cex.main = 1.5,
    cex.lab = 1.5,
    cex.axis = 1.5
  )
  
  return(sml)
}

pipeline <- function(cohort_num, txdb_regions, window_name) {
  if(cohort_num==1) {
    cohort_betas <- cohort1_betas
    cohort_covar <- cohort1_covars
    pltTitle <- "Cohort 1 (n=32)"
  } else if(cohort_num==2) {
    cohort_betas <- cohort2_betas
    cohort_covar <- cohort2_covars
    cohort_covar$Sample_ID <- cohort_covar$Tube.number
    pltTitle <- "Cohort 2 (n=58)"
  }
  
  TNBC_HRD <- cohort_covar$Sample_ID[cohort_covar$HRD == "Yes"]
  TNBC_CTRL <- cohort_covar$Sample_ID[cohort_covar$HRD == "No"]
  
  GR_HRD <- assemble_summ_beta_granges(cohort_betas[,TNBC_HRD], annot_genomation, rowMedians)
  GR_CTRL <- assemble_summ_beta_granges(cohort_betas[,TNBC_CTRL], annot_genomation, rowMedians)
  
  wrapper_genom(
    GR_HRD, 
    GR_CTRL, 
    txdb_regions, 
    window_name = window_name, 
    meth_var_name="median beta-value", pltTitle=pltTitle
  )
}

main <- function() {
  load_cohort1()
  load_cohort2()
  
  pdf("~/Downloads/HRD_txn_regulatory_methyl.pdf", width=8.27, height=8.27)
  for(cohort_num in c(1,2)) {
    ## Generic hg19:
    pipeline(cohort_num, TRANSCRIPT_REGIONS, "Distance (bp) to Transcription Start")

    ## HMEC ChIP-seq reference maps:
    pipeline(cohort_num, PEAKS_HMEC_H3K27ac, "Distance (bp) to H3K27ac")
    pipeline(cohort_num, PEAKS_HMEC_H3K36me3, "Distance (bp) to H3K36me3")
  }
  dev.off()
}

if(! interactive()) main()
