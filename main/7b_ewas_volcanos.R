# Volcano Plots of EWAS Results
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../utils.R")
source("../plotThemes.R")
library(ggrepel)
library(VennDiagram)

OUT_DIR <- "**************** MASKED ****************"

volcano_wrapper <- function(DMPs, xThreshUp=5.5, xThreshDown=-5.5, textHeight=3,
                            pThresh=0.05, esThresh=0) {
  #'@description Wrapper to draw volcano plot for differentially methylated CpGs
  DMPs$negLog10Pval <- -log10(DMPs$adj.P.Val)
  
  ## Label for direction of change:
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC > esThresh] <- "Positive"
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC < -esThresh] <- "Negative"
  DMPs$dir[is.na(DMPs$dir)] <- ""
  
  ## Add n each group to label:
  nPos <- sum(DMPs$dir=="Positive", na.rm=TRUE)
  nNeg <- sum(DMPs$dir=="Negative", na.rm=TRUE)
  
  ## Label selected genes:
  DMP.genes <- as.character(DMPs$UCSC_RefGene_Name)
  DMP.genes <- strsplit(DMP.genes, split=";")
  DMPs$Gene <- rep(NA, length(DMP.genes))
  for(k in 1:length(DMP.genes)) {
    DMP.genes[[k]] <- unique(DMP.genes[[k]])
    if(length(DMP.genes[[k]]) > 1) {
      DMPs$Gene[k] <- paste(DMP.genes[[k]], collapse=";")
    } else if(length(DMP.genes[[k]]) == 1) {
      DMPs$Gene[k] <- DMP.genes[[k]]
    } else if(length(DMP.genes[[k]]) == 0) {
      DMPs$Gene[k] <- NA
    }
  }
  DMPs$Label <- DMPs$Gene
  
  ## First pass for gene labels: remove genes taht didn't meet threshold:
  DMPs$Label[(DMPs$logFC <= xThreshUp & DMPs$logFC >= xThreshDown) | DMPs$adj.P.Val >= pThresh] <- NA #higher threshold for gene labeling
  
  ## (Optional) 2nd pass for gene labels: add back some:
  # cond_add_back <- DMPs$negLog10Pval > 1.25
  # DMPs$Label[cond_add_back] <- DMPs$Gene[cond_add_back]
  
  ## Auto determine text position & axis range/ticks:
  hMax <- max(DMPs$negLog10Pval)
  lMax <- min(DMPs$logFC)
  rMax <- max(DMPs$logFC)
  
  plt <- ggplot(DMPs, aes(x=logFC, y=negLog10Pval, color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_color_manual(values = c("gray","royalblue","red")) +
    scale_size_manual(values=c(1, 1.1, 1.1)) +
    scale_alpha_manual(values=c(0.9, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=-log10(pThresh), linetype="dashed") +
    labs(x="Log2 Fold Difference", y="-Log10 FDR") +
    annotate("text", rMax/2, hMax+0.5,  label=paste(nPos,"pos. association \n with HRD"), size=6, color="red") +
    annotate("text", lMax/2, hMax+0.5, label=paste(nNeg,"neg. association \n with HRD"), size=6, color="royalblue") +
    myVolcanoTheme
  
  if(esThresh == 0) {
    plt <- plt + geom_vline(xintercept=c(esThresh), linetype="dashed");
  } else {
    plt <- plt + geom_vline(xintercept=c(-esThresh, esThresh), linetype="dashed");
  }
  return(plt)
}

overlap_by_dir <- function(df1, dir1="UP", df2, dir2="UP", fname=NULL) {
  sele <- list(
    Cohort1 = df1$Name[df1$SignifByHRD %in% dir1],
    Cohort2 = df2$Name[df2$SignifByHRD %in% dir2]
  )
  venn.diagram(
    sele, 
    main = paste("Cohort1:", dir1, "Cohort2:", dir2),
    cex=1, cat.sex=2, cat.default.pos="text",
    margin = 0.05,
    disable.logging = TRUE,
    filename=paste0(OUT_DIR,fname)
  )
}

get_concordant <- function(df1, df2, dir_change) {
  df1 <- subset(df1, SignifByHRD==dir_change)
  df2 <- subset(df2, SignifByHRD==dir_change)
  
  res <- merge(df1, df2[ , COLS_EWAS], by="Name")
  colnames(res) <- gsub(".x", "_cohort1", colnames(res), fixed=TRUE)
  colnames(res) <- gsub(".y", "_cohort2", colnames(res), fixed=TRUE)
  res$log2FC_avg_cohorts <- (res$logFC_cohort1 + res$logFC_cohort2) / 2
  return(res[order(res$log2FC_avg_cohorts, decreasing=TRUE), ])
}

main <- function() {
  ewas_cohort1 <- read.csv(PATH_DMPS[["Cohort1"]], stringsAsFactors=FALSE)
  ewas_cohort2 <- read.csv(PATH_DMPS[["Cohort2"]], stringsAsFactors=FALSE)
  
  # -------------------- Volcano plots --------------------
  sum(ewas_cohort1$adj.P.Val < 0.05 & ewas_cohort1$logFC > 5.5)
  sum(ewas_cohort1$adj.P.Val < 0.05 & ewas_cohort1$logFC < -4.85)
  volcano_wrapper(ewas_cohort1, 5.5, -4.85, 2.5, esThresh=1)
  
  sum(ewas_cohort2$adj.P.Val < 0.05 & ewas_cohort2$logFC > 3.5)
  sum(ewas_cohort2$adj.P.Val < 0.05 & ewas_cohort2$logFC < -4.25)
  volcano_wrapper(ewas_cohort2, 3.5, -4.25, 5, esThresh=1)
  
  # ---------------------- Venn Diagrams ---------------------- 
  ewas_cohort1 <- ewas_df_cleanup(ewas_cohort1)
  ewas_cohort2 <- ewas_df_cleanup(ewas_cohort2)
  
  ## Concordant loci
  overlap_by_dir(ewas_cohort1, "UP", ewas_cohort2, "UP", "c1c2UP.png")
  overlap_by_dir(ewas_cohort1, "DOWN", ewas_cohort2, "DOWN", "c1c2DOWN.png")
  
  ## Discortant loci
  overlap_by_dir(ewas_cohort1, "UP", ewas_cohort2, "DOWN", "c1UP_c2DOWN.png")
  overlap_by_dir(ewas_cohort1, "DOWN", ewas_cohort2, "UP", "c1DOWN_c2UP.png")
  
  ## Concordant set: 639 loci
  ewas_concord <- rbind(
    get_concordant(ewas_cohort1, ewas_cohort2, "UP"),
    get_concordant(ewas_cohort1, ewas_cohort2, "DOWN")
  )
  write.csv(ewas_concord, PATH_DMPS[["Concordant"]], row.names=FALSE, quote=FALSE) #for heatmap/suppl.table
}

if(! interactive()) main()