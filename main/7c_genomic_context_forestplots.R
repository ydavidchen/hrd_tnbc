# Genomic Context Enrichment/Depletion Analysis of EWAS Loci

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggstance)
library(reshape2)
source("../utils.R")
source("../plotThemes.R")
source("../utils/ewas_helpers.R")

CONTEXTS <- c("Island","OpenSea", "Shore","Shelf")

genomic_enrich <- function(cohort_num, direction, contexts, confLev=NULL) {
  #'@description Wrapper to compute proportions in Significant vs. Background Sets
  #'#'@param direction: Direction of association
  #'@param confLev: Optional. If set, Fisher's test will be run
  
  dmps <- read.csv(PATH_DMPS[[paste0("Cohort_",cohort_num)]])
  dmps <- ewas_df_clean_up(dmps, 1.0, 0.05, FALSE)
  
  ## One-hot encode:
  if("Island" %in% contexts) dmps$Island <- dmps$Relation_to_Island == "Island"
  if("Shore" %in% contexts) dmps$Shore <- grepl("_Shore", dmps$Relation_to_Island)
  if("Shelf" %in% contexts) dmps$Shelf <- grepl("_Shelf", dmps$Relation_to_Island)
  if("OpenSea" %in% contexts) dmps$OpenSea <- dmps$Relation_to_Island == "OpenSea"

  ## Odds ratios for each genomic context:
  dmps_sig <- subset(dmps, SignifByHRD==direction)
  summFisher <- data.frame(Context=contexts, 
                           num_signif_yes=NA, num_signif_no=NA, prop_signif_yes=NA,
                           num_input_yes=NA, num_input_no=NA, prop_input_yes=NA, 
                           OR=NA, CI.lower=NA, CI.upper=NA, P=NA)
  for(gc in contexts) {
    contTab <- rbind(
      Signif = table(dmps_sig[gc]),
      Input = table(dmps[gc])
    )
    contTab <- t(contTab[ , c(2,1)])
    summFisher$num_signif_yes[summFisher$Context==gc]  <- contTab[1,1]
    summFisher$num_input_yes[summFisher$Context==gc]   <- contTab[1,2]
    summFisher$prop_signif_yes[summFisher$Context==gc] <- contTab[1,1] / sum(contTab[1:2,1])
    
    summFisher$num_signif_no[summFisher$Context==gc]  <- contTab[2,1]
    summFisher$num_input_no[summFisher$Context==gc]   <- contTab[2,2]
    summFisher$prop_input_yes[summFisher$Context==gc] <- contTab[1,2] / sum(contTab[1:2,2])
    
    if(!is.null(confLev)) {
      fT <- fisher.test(contTab, conf.level=confLev)
      summFisher$OR[summFisher$Context==gc]       <- fT$estimate
      summFisher$CI.lower[summFisher$Context==gc] <- fT$conf.int[[1]] 
      summFisher$CI.upper[summFisher$Context==gc] <- fT$conf.int[[2]]
      summFisher$P[summFisher$Context==gc] <- fT$p.value
    }
  }
  return(summFisher)
}

wrapper_forest <- function(cohort_num, res_fisher) {
  res_rf <- rbind(
    res_fisher[[paste0("c", cohort_num, "up")]],
    res_fisher[[paste0("c", cohort_num, "down")]]
  )
  res_rf$`Association with HRD` <- NA
  res_rf$`Association with HRD`[1:length(CONTEXTS)] <- "Pos"
  res_rf$`Association with HRD`[(length(CONTEXTS)+1):nrow(res_rf)] <- "Neg"
  
  res_rf$log2OR <- log2(res_rf$OR)
  res_rf$CI.lower <- log2(res_rf$CI.lower)
  res_rf$CI.upper <- log2(res_rf$CI.upper)
  
  ggplot(res_rf, aes(log2OR, Context, xmin=CI.lower, xmax=CI.upper)) +
    ggstance::geom_pointrangeh(aes(color=`Association with HRD`), size=1, position=position_dodge(width=0.5)) +
    geom_vline(xintercept=0, size=0.5, linetype="dashed") +
    labs(x="Log2 Odds Ratio (95% CI)") +
    scale_color_brewer(palette="Set1", direction=-1) +
    myForestTheme
}

main <- function() {
  RES_FISHER <- list(
    c1up = genomic_enrich(1, "UP", CONTEXTS, 0.95),
    c2up = genomic_enrich(2, "UP", CONTEXTS, 0.95),
    c1down = genomic_enrich(1, "DOWN", CONTEXTS, 0.95),
    c2down = genomic_enrich(2, "DOWN", CONTEXTS, 0.95)
  )
  print(RES_FISHER)
  
  wrapper_forest(1, RES_FISHER)
  wrapper_forest(2, RES_FISHER)
}

if(! interactive()) main()
