# Global Methylation Clusters w/ Recursively Partitioned Mixture Model

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")

SELE_UNIV <- read.table(PATH_UNIV)[,1] #25081 CpGs

getRPMMSampOrder <- function(rpmmClusters, Y_inv, sampleIdName) {
  #'@description Retrieves sample orders fo heat map visualization
  #'@param rpmmClusters data.frame with row.names = sample ID & 1 column named RPMM with cluster assignments
  #'@param Y_inv Input matrix for RPMM computation
  sampOrder <- c()
  for(r in names(table(rpmmClusters$RPMM))) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == r]
    clu <- t(Y_inv[rownames(Y_inv) %in% as.character(samps), ])
    s_i <- seriation::seriate(clu, margin=2)
    so_i <- seriation::get_order(s_i)
    sampOrder <- c(sampOrder, samps[so_i])
  }
  sampOrder <- data.frame(
    V1 = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  )
  colnames(sampOrder)[1] <- sampleIdName
  return(sampOrder)
}

run_custom_rpmm <- function(dnam, sele, maxlevel=2, sampleIdName="Sample_ID") {
  require(RPMM)
  Y_inv <- t(dnam[sele, ])
  sup_rpmm <- blcTree(Y_inv, verbose=1, maxlevel=maxlevel)
  print(sup_rpmm)
  plot(sup_rpmm)
  
  sup_class <- blcTreeLeafClasses(sup_rpmm)
  res_rpmm <- data.frame(table(RPMM=sup_class, Sample_ID=rownames(Y_inv)))
  res_rpmm <- subset(res_rpmm, Freq != 0)
  res_rpmm$Freq <- NULL
  
  rownames(res_rpmm) <- res_rpmm$Sample_ID #for sample order helper
  sampOrder <- getRPMMSampOrder(res_rpmm, Y_inv, sampleIdName)
  res_rpmm <- merge(res_rpmm, sampOrder, by=sampleIdName)
  return(res_rpmm)
}

main <- function() {
  load_cohort1()
  load_cohort2()
  
  c1_rpmm <- run_custom_rpmm(cohort1_betas, SELE_UNIV, "Sample_ID")
  c2_rpmm <- run_custom_rpmm(cohort2_betas, SELE_UNIV, "Tube.number")

  write.csv(c1_rpmm, paste0(DIR_RPMM,"Cohort1_RPMM.csv"), row.names=FALSE, quote=FALSE)
  write.csv(c2_rpmm, paste0(DIR_RPMM,"Cohort2_RPMM.csv"), row.names=FALSE, quote=FALSE)
}

if(! interactive()) main()
