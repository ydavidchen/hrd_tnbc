# Helper Functions
# Data Cleaning & Exploratory Analysis Procedures maybe masked

DIR_COHORTS <- "************ MASKED ************"
DIR_RPMM <- "************ MASKED ************"
COHORT1_PATHS <- "************ MASKED VECTOR ************"
COHORT2_PATHS <- "************ MASKED VECTOR ************"
DIR_SNPS <- "************ MASKED ************"
DIR_SYNAPSE <- "************ MASKED ************"

#------------------------------------- Data Loading Procedures -------------------------------------
load_new_rpmm <- function(cohort, dir=DIR_RPMM, suffix="_RPMM.csv") {
  rpmmDf <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  return(rpmmDf)
}

load_cohort1 <- function(dir=DIR_COHORTS) {
  DNAM_PATH <- " *********** MASKED ***********"
  INFO_PATH <- " *********** MASKED ***********"
  RAW_COVAR_PATH <- " *********** MASKED ***********"
  SELE_COLS <- " *********** MASKED ***********"
  
  ## Methylation data & sample info:
  cohort1_info <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  load(DNAM_PATH)
  
  cohort1_betas <- EPIC.betas[ , colnames(EPIC.betas) %in% cohort1_info$Sample_ID]
  
  ## Add new RPMM results calculated on Common CpG set
  c1_rpmm <- load_new_rpmm(1, paste0(dir,COHORT1_PATHS[4]))
  cohort1_info <- merge(cohort1_info, c1_rpmm, by="Sample_ID")
  
  ## Clean raw covariates:
  cohort1_covars <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  cohort1_covars <- cohort1_covars[ , SELE_COLS]
    
  ## Add in covariates & clean data:
  cohort1_covars <- merge(cohort1_info, cohort1_covars, by="Sample_ID")
  cohort1_covars <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  
  cohort1_betas <- cohort1_betas[ , match(cohort1_covars$Sample_ID,colnames(cohort1_betas))]
  stopifnot(identical(colnames(cohort1_betas), cohort1_covars$Sample_ID)) #checkpoint
  assign("cohort1_betas", cohort1_betas, 1) #1=.GlobalEnv
  assign("cohort1_covars", cohort1_covars, 1)
}

load_cohort2 <- function(dir=DIR_COHORTS, fname=COHORT2_PATHS[1], annotname=COHORT2_PATHS[2],
                         dir_rpmm=COHORT2_PATHS[3], path_addl=COHORT2_PATHS[4]) {
  ## Data loading & cleaning:
  load(paste0(dir,fname))
  cohort2_covars <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  
  ## Add latest RPMM results:
  c2_rpmm <- load_new_rpmm(2, dir=paste0(dir,dir_rpmm))
  cohort2_covars <- merge(cohort2_covars, c2_rpmm, by="Tube.number")
  
  ## Add in Age & Stage:
  cohort2_age_stage <- " *********** DATA CLEANING PROCEDURE MASKED ***********"
  cohort2_covars <- merge(cohort2_covars, cohort2_age_stage, by="Year_Patient_Block")
  
  ## Match & return:
  cohort2_betas <- cohort2_betas[ , match(cohort2_covars$Tube.number, colnames(cohort2_betas))]
  stopifnot(identical(colnames(cohort2_betas), cohort2_covars$Tube.number)) #checkpoint
  assign("cohort2_betas", cohort2_betas, 1)
  assign("cohort2_covars", cohort2_covars, 1)
}


loadSelectedCpGs <- function(dir, fname="***** MASKED ******") {
  cpgs <- read.table(paste0(dir,fname),stringsAsFactors=FALSE)[ , 1]
  return(cpgs)
}

loadEPICannotationFile <- function() {
  #'@description Loads MethylationEPIC 850K annotation as a data.frame
  require(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
  annot.850kb3 <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19))
  annot.850kb3$Methyl450_Loci <- annot.850kb3$Methyl450_Loci == "TRUE"
  annot.850kb3$Methyl27_Loci <- annot.850kb3$Methyl27_Loci == "TRUE"
  annot.850kb3$isEnhancer <- annot.850kb3$X450k_Enhancer=="TRUE" | annot.850kb3$Phantom4_Enhancers != "" | annot.850kb3$Phantom5_Enhancers != ""
  annot.850kb3$isPromoter <- grepl("TSS", annot.850kb3$UCSC_RefGene_Group)
  annot.850kb3[annot.850kb3==""] <- NA
  return(annot.850kb3)
}

#------------------------------------- Exploration & Statistical Methods -------------------------------------
custom_imputation <- function(mat, ...) {
  #'@description KNN imputation for DNA methylation using K=5 (Teschendorff 2016)
  if(mean(is.na(mat)) == 0) {
    warning("No missing values. Original matrix returned!");
    return(mat)
  }
  require(impute)
  temp <- impute.knn(mat, k=5, ...)
  return(temp$data)
}

display_contingency <- function(ctab, margin=2) {
  #'@description Assembles pretty contingency tables with row/column proportions
  #'@param ctab R matrix representing contingency table
  #'@param margin Whether % is to be calculated from row or column total
  ctab_prop <- round(100*prop.table(ctab,margin=margin), 2)
  stopifnot(identical(dimnames(ctab), dimnames(ctab_prop)))
  for(k in 1:length(ctab)) ctab[k] <- paste0(ctab[k], " (", ctab_prop[k], "%)")
  return(ctab)
}

run_2x2_assoc <- function(var1, var2, data, flipVar1=FALSE, flipVar2=FALSE, test="no") {
  contTab <- table(data[ , var1], data[ , var2], useNA="ifany")
  if(flipVar1) contTab <- contTab[c(2,1), ]
  if(flipVar2) contTab <- contTab[ , c(2,1)]
  print(display_contingency(contTab))
  if(test=="fisher") {
    fisher.test(contTab)
  } else if(test=="mcnemar") {
    mcnemar.test(contTab)
  }
}


fread_dnam <- function(csv_path, indColName="V1") {
  #'@description Fast-reads CSV files of GEO series matrix file beta-values re-saved as CSV files
  dnam <- data.table::fread(csv_path, data.table=FALSE, header=TRUE)
  rownames(dnam) <- dnam[ , colnames(dnam)==indColName]
  dnam <- dnam[ , -1]
  dnam <- data.matrix(dnam)
  return(dnam)
}

combine_dnam_matrices <- function(mat1, mat2) {
  #'@description Helper function to join/merge 2 DNA methylation matrices w/ rows=CpGs
  reformat_dnam_matrix <- function(mat) {
    mat <- as.data.frame(mat)
    mat$CpG <- rownames(mat)
    rownames(mat) <- NULL
    return(mat)
  }
  
  mega_data <- merge(
    reformat_dnam_matrix(mat1),
    reformat_dnam_matrix(mat2), 
    by = "CpG"
  )
  rownames(mega_data) <- mega_data$CpG
  mega_data$CpG <- NULL
  
  return(data.matrix(mega_data))
}

decideNumberOfMostVariable <- function(data, varThresh, plot=TRUE) {
  #'@description Subset a matrix of CpGs based on sample variance
  #'@param data Matrix of CpGs with rows=CpGs, columns=samples
  #'@param varThresh Variance threshold for color in red
  #'@param plot Should a ranked variance distribution be shown?
  require(matrixStats)
  
  vars <- matrixStats::rowVars(data)
  vars <- sort(vars, decreasing=TRUE)
  bool <- (vars >= varThresh)
  k <- sum(bool)
  
  if(plot) {
    plot(
      vars,
      col = ifelse(bool, "red", "black"),
      cex = 0.3,
      bty = "l",
      xlab = "CpGs",
      ylab = "Variance",
      main = "Inter-sample Variance Distribution"
    );
    abline(h=varThresh, lty=2)
    text(1.5e5, 0.10, paste(k, "CpGs"), col="red")
  }

  return(k)
}

selectMostVariableCpGs <- function(data, k) {
  #'@description Subset a matrix of CpGs based on sample variance
  #'@param data Matrix of CpGs with rows=CpGs, columns=samples
  #'@param k Number of most variable CpGs to select
  require(matrixStats)
  sele <- order(matrixStats::rowVars(data), decreasing=TRUE)[1:k]
  mat <- data[sele, ]
  return(mat)
}

createCpGTrackingBars <- function(stratifyPromoter=FALSE, noAsBlanks=FALSE) {
  #'@description Creates CpG annotation data.frame for tracking bars in a pheatmap
  annot.850kb3 <- loadEPICannotationFile()
  row_annot <- data.frame(
    row.names = annot.850kb3$Name, 
    Context = gsub("[N,S]_", "", annot.850kb3$Relation_to_Island),
    TSS200 = ifelse(grepl("TSS200", annot.850kb3$UCSC_RefGene_Group), "Yes", "No"),
    TSS1500 = ifelse(grepl("TSS1500", annot.850kb3$UCSC_RefGene_Group), "Yes", "No"),
    Enhancer = ifelse(annot.850kb3$isEnhancer, "Yes", "No")
  )
  row_annot$Promoter <- ifelse(row_annot$TSS200=="Yes" | row_annot$TSS1500=="Yes", "Yes", "No")
  if(! stratifyPromoter) row_annot$TSS200 <- row_annot$TSS1500 <- NULL
  if(noAsBlanks) row_annot[row_annot=="No"] <- NA
  return(row_annot)
}

