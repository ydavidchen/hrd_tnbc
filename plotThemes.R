# Reusable Objects & Constants for Data Visualization
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

library(ggplot2)

## Heat maps:
HEAT_COLS <- colorRampPalette(c("yellow","black","blue"))(1024)
HEAT_COLS_CNA <- colorRampPalette(c("blue","white","red"))(1024)

COLORS_BINARY <- c(Yes="black", No="lightgray")

ANNOT_COLORS <- list(
  `BRCA Mutation` = c(Yes="black", `No/unknown`="lightgray"),
  RPMM = c(A="black", B1="darkgray", B2="lightgray",
           C="black", D1="darkgray", D2="lightgray"),
  HRD = COLORS_BINARY,
  Mutation = c(BRCA1="red", BRCA2="purple"),
  methylBRCA1 = COLORS_BINARY,
  Promoter = COLORS_BINARY,
  Enhancer = COLORS_BINARY,
  Context = c(Island="black", Shore="dimgray", Shelf="gray60", OpenSea="gray90")
)

## ggplots
myScatterTheme <- theme_bw() + 
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        title = element_text(size=20,color="black"), 
        panel.spacing = unit(0.5, "lines"), panel.border = element_blank(), axis.line=element_line(color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black"))

myBarplotTheme <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") )

myBoxplotTheme <- theme_bw() +
  theme(axis.text.x=element_text(size=20,color="black"), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        strip.text.x=element_text(size=20,color="black",face="bold"), strip.background=element_rect(fill="gray95"),
        panel.border = element_blank(), axis.line=element_line(color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black"))

myBoxplotTheme2 <- myBoxplotTheme
myBoxplotTheme2$axis.text.x <- element_text(size=20, color="black")
myBoxplotTheme2$axis.title.x <- element_text(size=20, color="black")

myVolcanoTheme <- theme_classic() +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_text(size=21, color="black"),
        legend.title=element_blank(), legend.text=element_blank(), legend.position="none")

myForestTheme <- theme_bw() +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_blank(),
        legend.position="top", legend.title=element_text(size=16), legend.text=element_text(size=14))

mySurvTheme <- theme_classic() +
  theme(axis.text=element_text(size=20,color="black"), 
        axis.title=element_text(size=20,color="black"),
        title=element_text(size=20,color="black"),
        legend.title=element_text(size=16,color="black",face="bold"), legend.text=element_text(size=16,color="black",face="bold"), 
        legend.box.background=element_rect(colour="black",size=2),
        text=element_text(size=20),
        strip.text.x=element_text(size=20,colour="black",face="bold"))
