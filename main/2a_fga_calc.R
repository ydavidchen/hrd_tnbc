# Estimating Fraction of Genome Altered (FGA) from CopyNumber450KCancer Output
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
PATTERN <- "sample_|.seg"
SEG_DIR_C1 <- "******** MASKED ******** "
SEG_DIR_C2 <-  "******** MASKED ******** "
PEAKPATH_C1 <-  "******** MASKED ******** "
PEAKPATH_C2 <-  "******** MASKED ******** "
OUT_DIR <-  "******** MASKED ******** "

execute_fga <- function(seg_dir, peaks_path, thresh_segmean=0.1, thresh_nummark=20, thresh_pval=0.01) {
  ## Load Peak-shifts file:
  PEAK_SHIFTS <- read.csv(peaks_path, row.names="Sample", stringsAsFactors=FALSE)
  PEAK_SHIFTS <- PEAK_SHIFTS[ , "Shifting", drop=FALSE]
  
  ## Iterative peak shifting & filtering
  segs <- list.files(seg_dir)
  
  sample_stats <- data.frame(
    row.names = gsub(PATTERN, "", segs), #sample name
    nseg = rep(NA, length(segs)), 
    nseg_raw=NA, fga_raw=NA, 
    nseg_auto=NA, fga_auto=NA,
    nseg_auto_probes=NA, fga_auto_probes=NA,
    nseg_auto_probes_pval=NA, fga_auto_probes_pval=NA,
    stringsAsFactors = FALSE
  )
  
  for(seg in segs) {
    ## Load individual segmented CNA, and add a column for shifted mean:
    samplename <- gsub(PATTERN, "", seg)
    segment <- read.table(paste0(seg_dir, seg), header=TRUE, stringsAsFactors=FALSE) #reusable
    segment$mean.shifted <- segment$seg.mean + PEAK_SHIFTS[samplename, ]
    
    ## 1) Calculate FGA for all BEFORE autocorrections WITHOUT filters:
    GENERAL_NORMALIZER <- sum(segment$loc.end - segment$loc.start)
    segment_raw <- segment[abs(segment$seg.mean) >= thresh_segmean, ]
    fga_raw <- sum(segment_raw$loc.end - segment_raw$loc.start) / GENERAL_NORMALIZER
    
    ## 2) Calculate FGA AFTER autocorrection WITHOUT filters: 
    segment_auto <- segment[abs(segment$mean.shifted) >= thresh_segmean, ]
    fga_auto <- sum(segment_auto$loc.end - segment_auto$loc.start) / GENERAL_NORMALIZER
    
    ## Compute a Common Denominator based on number of probes:
    segment_probes <- segment[segment$num.mark >= thresh_nummark, ] #reusable
    PROBE_BASED_NORMALIZER <- sum(segment_probes$loc.end - segment_probes$loc.start) #reusable
    
    ## 3) Filter by autocorrected segment mean AND number of probes: 
    segment_auto_probes <- segment_probes[abs(segment_probes$mean.shifted) >= thresh_segmean, ]
    fga_auto_probes <- sum(segment_auto_probes$loc.end - segment_auto_probes$loc.start) / PROBE_BASED_NORMALIZER
    
    ## 4) Filter by autocorrected segment mean AND number of probes AND P-value:
    segment_auto_probes_pval <- segment_auto_probes[segment_auto_probes$pval <= thresh_pval, ]
    fga_auto_probes_pval <- sum(segment_auto_probes_pval$loc.end - segment_auto_probes_pval$loc.start, na.rm=TRUE) / PROBE_BASED_NORMALIZER
    
    ## Iteratively combine into data.frames by sample:
    sample_stats[samplename, ] <- c(
      nrow(segment), 
      nrow(segment_raw), fga_raw,
      nrow(segment_auto), fga_auto,
      nrow(segment_auto_probes), fga_auto_probes,
      sum(! is.na(segment_auto_probes_pval$pval)), fga_auto_probes_pval
    )
  }
  return(sample_stats)
}

fga_cohort1 <- execute_fga(SEG_DIR_C1, PEAKPATH_C1)
fga_cohort2 <- execute_fga(SEG_DIR_C2, PEAKPATH_C2)

# write.csv(fga_cohort1, file=paste0(OUT_DIR,"cohort1_fga_stats.csv"), row.names=TRUE, quote=FALSE)
# write.csv(fga_cohort2, file=paste0(OUT_DIR,"cohort2_fga_stats.csv"), row.names=TRUE, quote=FALSE)
