## Table One Summary: Overall & Within BRCA1-altered

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
library(tableone)

## Load Data:
load_cohort1(); load_cohort2()
rm(cohort1_betas, cohort2_betas)

# --------------------------------- TNBC Cohort 1 ---------------------------------
## Overall:
print(CreateTableOne(
  c("AgeAtDx", "binaryStage", "BRCA.Status", "BRCA1_Hypermeth"), 
  strata = "MLPA",
  data = cohort1_covars,
  includeNA = TRUE,
  test = FALSE
), showAllLevels=TRUE)


## Within BRCA1-altered (Reviewer Request):
cohort1_covars <- subset(cohort1_covars, HRD=="Yes")
cohort1_covars <- subset(cohort1_covars, BRCA.Status=="BRCA1" | BRCA1_Hypermeth=="Yes")

print(CreateTableOne(
  c("AgeAtDx", "binaryStage", "MLPA.mean", "RPMM"), 
  strata = "BRCA1_Hypermeth",
  data = cohort1_covars, 
  includeNA = TRUE,
  test = FALSE
), showAllLevels=TRUE)

# --------------------------------- TNBC Cohort 2 ---------------------------------
## Overall:
print(CreateTableOne(
  c("AgeAtDx", "StageBinary", "BRCA.Status", "isBRCA1Meth"), 
  strata = "MLPA", 
  data = cohort2_covars, 
  includeNA = TRUE,
  test = FALSE
), showAllLevels=TRUE)


## Within BRCA1-altered (Reviewer Request):
cohort2_covars <- subset(cohort2_covars, HRD=="Yes")
cohort2_covars <- subset(cohort2_covars, BRCA.Status=="BRCA1" | isBRCA1Meth)
print(CreateTableOne(
  c("AgeAtDx", "StageBinary", "MLPA.mean", "RPMM"), 
  strata = "isBRCA1Meth", 
  data = cohort2_covars, 
  includeNA = TRUE,
  test = FALSE
), showAllLevels=TRUE)
