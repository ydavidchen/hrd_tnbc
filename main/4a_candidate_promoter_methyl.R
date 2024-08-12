# Hypermethylation Calls of High-profile Candidate Genes

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")

## Constants: see original publications Tutt 2018, Sharma 2014, Polak 2017
LOCI_TUTT <- c("cg24806953", "cg16630982", "cg16963062", "cg15419295", 
               "cg21253966", "cg20187250", "cg04110421","cg17301289", "cg04658354")

LOCI_SHARMA <- c("cg19531713", "cg08993267", "cg04658354") #cg19088651 didn't pass all QCs

LOCI_RAD51C <- c("cg27221688", "cg02118635", "cg24099023", "cg10487724", "cg07329131", "cg02110529", 
                 "cg14066645", "cg23436779", "cg16398362", "cg08140274", "cg06855092", "cg19208681", 
                 "cg25339112", "cg05214530", "cg27281549", "cg12544923")


THRESH_SHARMA <- THRESH_POLAK <- 0.2

THRESH_TUTT <- 0.1

KEY <- "Sample_ID"

## Reusable wrappers:
define_hyperme <- function(betas, loci, definition, key) {
  df <- data.frame(t(betas[rownames(betas) %in% loci, , drop=FALSE]))
  colnames(df) <- paste(definition, colnames(df), sep="_")
  
  if(definition == "Sharma") {
    df$SharmaCpGAvg <- rowMeans(df)
    df$SharmaBRCA1meth <- df$SharmaCpGAvg >= THRESH_SHARMA
  } else if(definition == "Tutt") {
    df$TuttNumOverThresh <- rowSums(df >= THRESH_TUTT)
    df$TuttBRCA1meth <- df$TuttNumOverThresh == length(loci)
  } else if(definition %in% c("RAD51C","Polak")) {
    df$RAD51COverThresh <- rowSums(df >= THRESH_POLAK) 
    df$RAD51Cmeth <- df$RAD51COverThresh == length(loci)
  }
  
  df[key] <- rownames(df)
  rownames(df) <- NULL
  return(df)
}

calc_dnam_statuses <- function(betas, key=KEY) {
  mTutt <- define_hyperme(betas, LOCI_TUTT, "Tutt", key)
  mSharma <- define_hyperme(betas, LOCI_SHARMA, "Sharma", key)
  mRad51c <- define_hyperme(betas, LOCI_RAD51C, "RAD51C", key)
  
  df <- mTutt
  df <- merge(df, mSharma, by=key)
  df <- merge(df, mRad51c, by=key)
  return(df)
}

# ------------------------ Impose Definitions ------------------------
load_cohort1()
load_cohort2()

c1m <- calc_dnam_statuses(cohort1_betas, "Sample_ID")
c2m <- calc_dnam_statuses(cohort2_betas, "Tube.number")

# write.csv(c1m, paste0(DIR_COHORTS,"BRCA1_methylation/Cohort1_methylated_samples.csv"), row.names=FALSE, quote=FALSE)
# write.csv(c2m, paste0(DIR_COHORTS,"BRCA1_methylation/Cohort2_methylated_samples.csv"), row.names=FALSE, quote=FALSE)

# ------------------------ Statistical Tests ------------------------
## BRCA1 methylation Assessment:
run_2x2_assoc("TuttBRCA1meth", "SharmaBRCA1meth", c1m, TRUE, TRUE, "none")
run_2x2_assoc("TuttBRCA1meth", "SharmaBRCA1meth", c2m, TRUE, TRUE, "none")

# Associations:
BRCAs <- c("BRCA1","BRCA2")

cohort1_covars <- merge(cohort1_covars, c1m, by="Sample_ID")
cohort1_covars$BRCA1Mut <- cohort1_covars$BRCA.Status == "BRCA1"
cohort1_covars$Altered <- cohort1_covars$BRCA.Status %in% BRCAs | cohort1_covars$TuttBRCA1meth

run_2x2_assoc("Altered", "HRD", subset(cohort1_covars, !is.na(BRCA.Status)), TRUE, FALSE, "none")
run_2x2_assoc("BRCA1Mut", "TuttBRCA1meth", subset(cohort1_covars, BRCA1Mut|TuttBRCA1meth), TRUE, FALSE, "none")

cohort2_covars <- merge(cohort2_covars, c2m, by="Tube.number")
cohort2_covars$BRCA1Mut <- cohort2_covars$BRCA.Status == "BRCA1"
cohort2_covars$Altered <- cohort2_covars$BRCA.Status %in% BRCAs | cohort2_covars$TuttBRCA1meth
run_2x2_assoc("Altered", "HRD", subset(cohort2_covars, !is.na(BRCA.Status)), TRUE, FALSE, "none")
run_2x2_assoc("BRCA1Mut", "TuttBRCA1meth", subset(cohort2_covars, BRCA1Mut|TuttBRCA1meth), TRUE, FALSE, "none")
