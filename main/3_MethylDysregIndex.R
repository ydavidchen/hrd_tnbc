# Methylation Dysregulation Index (MDI)
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
library(matrixStats)
library(reshape2)
source("../utils.R")

source("../plotThemes.R")
myBoxplotTheme$axis.title.x <- element_text(size=20,color="black")
myScatterTheme$panel.spacing.x <- ggplot2::unit(2, "lines")
KEY <- "Sample_ID"

## TNBCs & controls methylation data & covariates:
load_cohort1()
load_cohort2()
load(paste0(DIR_MDI,"Oltra2018HealthySubjects.RData"))

mdi_wrapper <- function(betas, healthyBetas=OltraHealthyBetas, universe=NULL) {
  #'@references O'Sullivan et al. Epigenet; Salas et al. Epigenet
  if(is.null(universe)) universe <- intersect(rownames(betas), rownames(healthyBetas))
  
  healthyBetas <- subset(healthyBetas, rownames(healthyBetas) %in% universe)
  betas <- subset(betas, rownames(betas) %in% universe)
  stopifnot(identical(rownames(betas), rownames(healthyBetas))) #checkpoint
  
  print(paste(length(universe), "CpGs used for computing MDI..."));
  absDeviats <- 100 * abs(betas - rowMedians(healthyBetas));
  mdi <- data.frame(
    Sample_ID = colnames(betas),
    MDI = colMeans2(absDeviats)
  )
  return(mdi)
}

stats_wrapper <- function(data, corMethod="pearson") {
  #'@description Wrapper to test MDI between HRD strata & correlation with MLPA score
  print(t.test(MDI ~ HRD, data))
  print("---------------------------------")
  print(cor.test(data$MDI, data$MLPA.mean, method=corMethod))
}

## Estimation & statistical tests:
mdi_c1 <- mdi_wrapper(cohort1_betas, OltraHealthyBetas, NULL) #823086 CpGs used for computing MDI...
mdi_c2 <- mdi_wrapper(cohort2_betas, OltraHealthyBetas, NULL) #733664 CpGs used for computing MDI...

# write.csv(mdi_c1, paste0(DIR_MDI,"200203_mdi_cohort1.csv"), row.names=FALSE, quote=FALSE)
# write.csv(mdi_c2, paste0(DIR_MDI,"200203_mdi_cohort2.csv"), row.names=FALSE, quote=FALSE)

## Join in clinical features:: 
mdi_c1 <- merge(mdi_c1, cohort1_covars[ , c(KEY,"MLPA.mean","HRD")], by=KEY)
mdi_c2 <- merge(mdi_c2, cohort2_covars[ , c("Tube.number","MLPA.mean","HRD")], by.x=KEY, by.y="Tube.number")

## Statistical tests:
stats_wrapper(mdi_c1)
stats_wrapper(mdi_c2)

## Publication-quality visualizations: 
mdi_cohorts <- rbind(
  cbind(Cohort="Cohort 1", mdi_c1),
  cbind(Cohort="Cohort 2", mdi_c2)
)

# png("~/Downloads/mdi_by_hrd_strata.png", width=5, height=10, units="in", res=300)
ggplot(mdi_cohorts, aes(HRD, MDI)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA, outlier.size=NA) +
  geom_jitter(width=0.25) +
  ylab("% Methylation Dysregulation") +
  annotate("text", 1.5, 18.5, label="*", size=12, fontface=2) +
  facet_wrap(~ Cohort) +
  myBoxplotTheme
# dev.off()

# png("~/Downloads/mdi_vs_hrdprob.png", width=15, height=7.5, units="in", res=300)
ggplot(mdi_cohorts, aes(MLPA.mean, MDI)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~ Cohort) +
  labs(x="HRD Probability", y="% Methylation Dysregulation Index") +
  geom_text(data=data.frame(x=0.35,y=16,Cohort="Cohort 1",lab="PCC (95% CI) = 0.30 (-0.06,0.58) \n P = 0.10"), aes(x=x,y=y,label=lab), size=5) +
  geom_text(data=data.frame(x=0.35,y=16,Cohort="Cohort 2",lab="PCC (95% CI) = 0.25 (-0.01,0.48) \n P = 0.06"), aes(x=x,y=y,label=lab), size=5) +
  myScatterTheme
# dev.off()
