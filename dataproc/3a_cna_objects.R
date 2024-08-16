# MethylationEPIC-basd CN Extraction & Initial Segmentation
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
library(minfi)
library(conumee)

CHROMS <- c(paste0("chr",1:22))
SEG_DIR_C1 <- "************** MASKED **************"
SEG_DIR_C2 <- "************** MASKED **************"

get_whole_genome_detailed_regions <- function(chroms=CHROMS) {
  #'@description GenomicRanges. Need to be formatted exactly the same as `data(detailed_regions)`
  #'@references https://genomicsclass.github.io/book/pages/bioc1_annoOverview.html
  require(GenomicRanges)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  
  ## "We can use genes() to get the addresses of genes using Entrez Gene IDs" 
  detailed_regions <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene, columns="gene_id")
  names(detailed_regions@elementMetadata) <- "name" #conumee requirement
  values(detailed_regions)$thick <- IRanges(
    start(detailed_regions) - (500000 - round(width(detailed_regions)/2)), 
    end(detailed_regions) + (500000 - round(width(detailed_regions)/2))
  );  #conumee requirement
  detailed_regions <- detailed_regions[detailed_regions@seqnames %in% chroms]
  return(detailed_regions)
}

assemble_annot_for_conumee <- function(detailed_regions, exclude_regions, chrXY=FALSE, array_type="EPIC") {
  #'@description: Creates hg19 annotation for conumee
  require(conumee);
  cAnnot <- CNV.create_anno(
    detail_regions = detailed_regions,
    exclude_regions = exclude_regions,
    array_type = array_type, 
    chrXY = chrXY
  )
  return(cAnnot)
}

cna_inference <- function(loadedChannel, ctrlChannel, conumeeAnno, pathForSegs=NULL) {
  #'@description Wrapper to infer Segmented CN iteratively
  #'@return List of list with gene level & segmented CNA sub-lists
  #'@details Saves un-processed segments to `pathForSegs`
  require(minfi)
  geneCNAList <- segmentList <- list()
  for(k in 1:length(names(loadedChannel))){
    sampName <- names(loadedChannel)[k]
    x <- CNV.fit(loadedChannel[k], ctrlChannel, anno=conumeeAnno)
    x <- CNV.bin(x)
    x <- CNV.detail(x)
    x <- CNV.segment(x)
    
    ## Store objects:
    geneCNAList[[sampName]] <- CNV.write(x, what="detail")
    segmentList[[sampName]] <- CNV.write(x, what="segments")
    if(! is.null(pathForSegs)) CNV.write(x, what="segments", file=paste0(pathForSegs,"sample_",sampName,".seg"))
  }
  return(list(geneCNAList=geneCNAList, segmentList=segmentList))
}

load_idat_for_cna_wrapper <- function(idat_dir, cohort_name) {
  #'@description Wrapper to load IDATs
  targets <- read.metharray.sheet(idat_dir)
  
  ## Customize existing tumor annotations based on cohort:
  if(cohort_name %in% c("2","cohort2")) {
    targets$Sample_ID <- targets$Tube.number
  } else if(cohort_name %in% c("control","healthy")) {
    targets$Sample_ID <- targets$Accession
  } else {
    error("Please name a valid cohort!")
  }
  methylObj <- read.metharray.exp(idat_dir, targets=targets) #RGChannelSet
  methylObj <- preprocessIllumina(methylObj) #Mset
  return(conumee::CNV.load(methylObj, names=pData(methylObj)$Sample_ID))
}

# ------------------------------ Build Shared Objects ------------------------------
## Assemble genome-scale annotations:
DETAILED_REGIONS <- get_whole_genome_detailed_regions()
data(exclude_regions) #conumee
CONUMEE_ANNOT <- assemble_annot_for_conumee(DETAILED_REGIONS, exclude_regions)

## Load IDAT files as channels:
control_ch <- load_idat_for_cna_wrapper(DIR_COHORTS[["Oltra"]], "control")
cohort1_ch <- load_idat_for_cna_wrapper(DIR_IDATE[["Cohort1"]], "cohort1")
cohort2_ch <- load_idat_for_cna_wrapper(DIR_IDATE[["Cohort2"]], "cohort2")

## Perform Probe overlap across cohorts, control & annotation:
common_probes <- Reduce(
  intersect,
  list(rownames(cohort1_ch@intensity), rownames(cohort2_ch@intensity), 
       rownames(control_ch@intensity), names(CONUMEE_ANNOT@probes))
)
length(common_probes)

CONUMEE_ANNOT@probes <- CONUMEE_ANNOT@probes[names(CONUMEE_ANNOT@probes) %in% common_probes, ]
control_ch@intensity <- control_ch@intensity[rownames(control_ch@intensity) %in% common_probes, ]
cohort1_ch@intensity <- cohort1_ch@intensity[rownames(cohort1_ch@intensity) %in% common_probes, ]
cohort2_ch@intensity <- cohort2_ch@intensity[rownames(cohort2_ch@intensity) %in% common_probes, ]

## Execute Inference: 
cna_obj_cohort1 <- cna_inference(cohort1_ch, control_ch, CONUMEE_ANNOT, SEG_DIR_C1)
cna_obj_cohort2 <- cna_inference(cohort2_ch, control_ch, CONUMEE_ANNOT, SEG_DIR_C2)
# save(list=c("cna_obj_cohort1","cna_obj_cohort2"), file=paste0(DIR_COHORTS,"CNA_raw_objects.RData"), compress=TRUE)
