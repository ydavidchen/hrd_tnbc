# Empirical Search for Common Universe of 25081 CpGs w/ Visualizations

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R") #imports helpers + constants
source("../plotThemes.R") #constants for graphics
myScatterTheme$axis.text.x <- element_text(size=13, angle=45)

find_optim_num <- function(vec_k) {
  #'@description Wrapper to iteratively examine 2-cohort intersection
  vec_count <- c()
  for(k in vec_k) {
    c1_univ <- rownames(selectMostVariableCpGs(cohort1_betas, k))
    c2_univ <- rownames(selectMostVariableCpGs(cohort2_betas, k))
    vec_count <- c(vec_count, length(intersect(c1_univ, c2_univ)))
  }
  
  plt <- ggplot() +
    geom_point(aes(vec_k, vec_count), size=2) +
    scale_x_continuous(breaks=vec_k) +
    labs(x="Number of most variable CpGs in each cohort", 
         y="Number of shared CpGs between cohorts") +
    geom_hline(yintercept=2.5e4, linetype="dashed", color="red") +
    myScatterTheme;
  
  print(plt);
  return(data.frame(n_queried=vec_k, n_shared=vec_count))
}

## Load process DNAm data:
load_cohort1()
load_cohort2()

## Export shared CpGs:
find_optim_num(seq(5000, 1e5, by=5000))
find_optim_num(seq(59850, 59950, by=5)) #zoom-in

## Genomic context distribution: 60K CpGs in each cohort:
K <- 6e4
seleUniv <- intersect(
  rownames(selectMostVariableCpGs(cohort1_betas, K)),
  rownames(selectMostVariableCpGs(cohort2_betas, K))
)

write.table(seleUniv, file=OUT_LIST_PATH, row.names=FALSE, quote=FALSE, col.names=FALSE)

