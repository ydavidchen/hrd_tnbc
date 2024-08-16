# Further Process & Normalize Segmented CNA from Methylation Arrays
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
library(CopyNumber450kCancer)

SEG_DIR_C1 <- "************** MASKED **************"
SEG_DIR_C2 <- "************** MASKED **************"
OUT_DIR_C1 <- "************** MASKED **************"
OUT_DIR_C2 <-  "************** MASKED **************"

SEGFILE_COL_CLASSES <- c("character","character","integer", "integer", "integer","numeric","numeric","numeric","numeric");

DATA_DICT <- data.frame(
  original = c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"),
  new = c("Sample","Chromosome","bp.Start","bp.End","Num.of.Markers","Mean")
)

combine_segments <- function(segment_dir) {
  #'@description Reads & combines segmented CNA files into a data.frame
  ## Iteratively load segmented CNA:
  files <- list.files(segment_dir)
  mergedSegments <- NULL
  for(f in files){
    file_loaded <- read.table(paste0(segment_dir, f), header=TRUE, stringsAsFactors=FALSE, colClasses=SEGFILE_COL_CLASSES)
    mergedSegments <- rbind(mergedSegments, file_loaded)
  }
  ## Swap column names based on data dictionary:
  mergedSegments <- mergedSegments[ , colnames(mergedSegments) %in% DATA_DICT$original]
  stopifnot(identical(DATA_DICT$original, colnames(mergedSegments)))
  colnames(mergedSegments) <- DATA_DICT$new
  return(mergedSegments)
}

peak_correction <- function(mergedSegments, out_dir) {
  #'@description (Void) Method to correct peaks
  ## Generate sample sheet for segment files:
  curr_dir <- getwd()
  
  sampleFileEpic <- data.frame(
    Number = 1:length(unique(mergedSegments$Sample)),
    Sample = unique(mergedSegments$Sample),
    Comment = "Case"
  )
  
  ## Generate a CopyNumber450kCancer object using segment files:
  setwd(out_dir) 
  object <- ReadData(
    regions_file = mergedSegments, 
    Sample_list = sampleFileEpic
  )
  object <- AutoCorrectPeak(object)
  PlotCNV(object, comment=FALSE)
  
  setwd(curr_dir) #reset
}

wrapper <- function(segment_dir, out_dir) {
  mergedSegments <- combine_segments(segment_dir);
  peak_correction(mergedSegments, out_dir);
}

wrapper(SEG_DIR_C1, OUT_DIR_C1)
wrapper(SEG_DIR_C2, OUT_DIR_C2)
