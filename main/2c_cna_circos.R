# Segmented CNA Visualization with Circos Plots
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(copynumber)
library(pheatmap)
source("../utils.R")
source("../plotThemes.R");

CHROMS <- paste0("chr", 1:22)
CN_COLS <- c(loss="red", gain="darkgreen")
THRESH_GAIN <- 0.10

SEG_DIR_C1 <- "********** MASKED **********"
SEG_DIR_C2 <- "********** MASKED **********"


assemble_long_seg_cn <- function(segment_dir) {
  #'@description Reads & assembles segmented CNA for `copynumber` package visualization
  seg_files <- list.files(segment_dir)
  df <- NULL
  for(fname in seg_files) {
    mat_k <- read.table(paste0(segment_dir, fname), header=TRUE)
    df <- rbind(df, mat_k)
  }
  df <- as.data.frame(df)
  df <- data.frame(
    sampleID = as.character(df$ID),
    chrom = as.factor(as.numeric(gsub("chr", "", df$chrom))),
    arm = NA,
    start.pos = df$loc.start,
    end.pos = df$loc.end,
    n.probes = df$num.mark,
    mean = df$seg.mean
  )
  return(df)
}

assemble_fga_matrix <- function(segment_dir, agg_fun=mean) {
  #'@description Aggregate segmented CNA (stored in individual text files per sample) and bind into matrix
  seg_files <- list.files(segment_dir)
  
  fga_mat <- data.frame(chrom=CHROMS)
  for(fname in seg_files) {
    mat_k <- read.table(paste0(segment_dir, fname), header=TRUE)
    mat_k <- aggregate(seg.mean ~ chrom, data=mat_k, FUN=agg_fun)
    colnames(mat_k)[colnames(mat_k)=="seg.mean"] <- fname
    fga_mat <- merge(fga_mat, mat_k, by="chrom")
  }
  rownames(fga_mat) <- fga_mat$chrom
  fga_mat$chrom <- NULL
  
  fga_mat <- fga_mat[match(CHROMS, rownames(fga_mat)), ]
  colnames(fga_mat) <- gsub("sample_", "", colnames(fga_mat))
  colnames(fga_mat) <- gsub(".seg", "", colnames(fga_mat))
  return(t(fga_mat))
}

## Load cohort annotations:
load_cohort1()
load_cohort2()
rm(cohort1_betas, cohort2_betas)

HRD_TNBCS <- c(
  cohort1_covars$Sample_ID[cohort1_covars$HRD == "Yes"], 
  cohort2_covars$Tube.number[cohort2_covars$HRD == "Yes"]
)

NONHRD_TNBCS <- c(
  cohort1_covars$Sample_ID[cohort1_covars$HRD == "No"], 
  cohort2_covars$Tube.number[cohort2_covars$HRD == "No"]
)

## Cohort segments:
cohort1_segs <- assemble_long_seg_cn(SEG_DIR_C1)
cohort2_segs <- assemble_long_seg_cn(SEG_DIR_C2)
cohorts_segs <- rbind(cohort1_segs, cohort2_segs)

## Circos Plots (Fig.1A): 
plotCircle(subset(cohorts_segs, sampleID %in% HRD_TNBCS), freq.colors=CN_COLS, thres.gain=THRESH_GAIN)
plotCircle(subset(cohorts_segs, sampleID %in% NONHRD_TNBCS), freq.colors=CN_COLS, thres.gain=THRESH_GAIN)

## Linear Gain/Loss Frequencies (Fig.S1): 
plotFreq(
  segments = subset(cohort1_segs, sampleID %in% HRD_TNBCS), 
  thres.gain = THRESH_GAIN, 
  col.gain = CN_COLS["gain"],
  col.loss = CN_COLS["loss"],
  main = "Cohort 1: HRD"
)

plotFreq(
  segments = subset(cohort1_segs, sampleID %in% NONHRD_TNBCS), 
  thres.gain = THRESH_GAIN,
  col.gain = CN_COLS["gain"],
  col.loss = CN_COLS["loss"],
  main = "Cohort 1: Non-HRD"
)

plotFreq(
  segments = subset(cohort2_segs, sampleID %in% HRD_TNBCS), 
  thres.gain = THRESH_GAIN, 
  col.gain = CN_COLS["gain"],
  col.loss = CN_COLS["loss"],
  main = "Cohort 2 (FFPE): HRD"
)

plotFreq(
  segments = subset(cohort2_segs, sampleID %in% NONHRD_TNBCS), 
  thres.gain = THRESH_GAIN, 
  col.gain = CN_COLS["gain"],
  col.loss = CN_COLS["loss"],
  main = "Cohort 2 (FFPE): Non-HRD"
)
