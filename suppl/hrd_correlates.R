# Molecular & Clinical Correlates of High vs. Low HRD Groups in TCGA

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(reshape2)
library(ggpubr)
library(survival)
library(survminer)

source("../utils.R")
source("../plotThemes.R")
myBoxplotTheme$axis.title.x <- element_text(size=20,color="black")
myBarplotTheme$axis.text.x <- element_blank()

PATH_BRCA <-  "**************** MASKED ****************" #TCGA TNBC BRCA1/2 mutations
PATH_SBS <- "**************** MASKED ****************" #curated from Rosenthal 2016
PATH_MUTS <- "**************** MASKED ****************" #curated from Kandoth 2013
PATH_SURV <-  "**************** MASKED ****************" #curated TCGA survival data
PATH_PSVR <- "**************** MASKED ****************" #sklearn output

ADMIN_CENSOR_Y <- 3
KEY <- "Tumor_Sample_Barcode"

## Supervised ML results:
FINAL_PREDS_TCGA <- read.csv(PATH_PSVR, stringsAsFactors=FALSE)
FINAL_PREDS_TCGA[KEY] <- FINAL_PREDS_TCGA$bcr_patient_barcode
FINAL_PREDS_TCGA$PSVR_strata <- ifelse(FINAL_PREDS_TCGA$PSVR > median(FINAL_PREDS_TCGA$PSVR), "High", "Low")

# ----------------- Load & clean published HRD-related metrics for TCGA -----------------
TCGABRCA_MSIG3 <- read.csv(PATH_SBS)

KANDOTH_MUTS <- gdata::read.xls(PATH_MUTS, stringsAsFactors=FALSE)
KANDOTH_MUTS$Tumor_Sample_Barcode <- substr(KANDOTH_MUTS$TCGA.ID, 1, 12)
KANDOTH_MUTS$Tumor_Sample_Barcode <- gsub(".", "-", KANDOTH_MUTS$Tumor_Sample_Barcode, fixed=TRUE)

TCGA_SURV <- read.csv(PATH_SURV, stringsAsFactors=FALSE)
colnames(TCGA_SURV)[1] <- KEY

SZTUPINSKI_BRCA <- read.csv(PATH_BRCA, stringsAsFactors=FALSE)
colnames(SZTUPINSKI_BRCA)[1] <- KEY
SZTUPINSKI_BRCA$HRD <- SZTUPINSKI_BRCA$HRDsum_WXS > 42
SZTUPINSKI_BRCA$isBRCAmut <- apply(SZTUPINSKI_BRCA[,c("BRCA1germline","BRCA2germline","BRCA1somatic","BRCA2somatic")], 1, FUN=function(vec) any(vec==1, na.rm=TRUE));
SZTUPINSKI_BRCA$isBRCAmut <- ifelse(SZTUPINSKI_BRCA$isBRCAmut, "BRCA1/2 mutated", "BRCA1/2 intact");

# ----------------- Association Analysis with Known HRD Correlates (Suppl Fig.) -----------------
## PSVR score profile:
plt_profile <- merge(FINAL_PREDS_TCGA, SZTUPINSKI_BRCA, by=KEY, all.x=TRUE)
plt_profile$isBRCAmut[is.na(plt_profile$isBRCAmut)] <- "(unknown)"
ggplot(plt_profile, aes(reorder(Tumor_Sample_Barcode,PSVR), PSVR, color=isBRCAmut, fill=isBRCAmut)) +
  geom_bar(stat="identity") +
  labs(x="113 TCGA TNBCs", y="PSVR output") +
  geom_vline(xintercept=57.5, linetype="dashed") +
  scale_color_manual(values=c("gray","black","red")) +
  scale_fill_manual(values=c("gray","black","red")) +
  annotate("text", 28, 1.05, label="Low HRD stratum", size=8) +
  annotate("text", 88, 1.05, label="High HRD stratum", size=8) +
  myBarplotTheme

## Compare/correlate w/ Telli HRD score:
t.test(HRDsum_WXS ~ PSVR_strata, data=plt_profile);
ggplot(plt_profile, aes(PSVR_strata, HRDsum_WXS)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA) +
  geom_jitter(width=0.25) +
  labs(x="Predicted HRD Stratum (PSVR)", y="Telli HRD Score") +
  annotate("text", 1.5, 100, label="** P = 0.004", size=8) +
  myBoxplotTheme

cor.test(plt_profile$HRDsum_WXS, plt_profile$PSVR);
ggplot(plt_profile, aes(HRDsum_WXS, PSVR)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Telli HRD Score", y="PSVR output") +
  scale_y_continuous(limits=c(0,1.1)) +
  annotate("text", 80, 0.25, label="PCC (95% CI) = 0.32 (0.11, 0.50) \n P = 0.003 **", size=5) +
  myScatterTheme

## Mutational Signature 3:
plt_cor_sig3 <- merge(FINAL_PREDS_TCGA, TCGABRCA_MSIG3, by=KEY)

t.test(100*Contribution ~ PSVR_strata, data=plt_cor_sig3)
ggplot(plt_cor_sig3, aes(PSVR_strata, 100*Contribution)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA) +
  geom_jitter(width=0.25) +
  annotate("text", 1.5, 70, label="* P = 0.024", size=8) +
  labs(x="Predicted HRD Stratum (PSVR)", y="Signature 3 contribution (%)") +
  myBoxplotTheme

cor.test(plt_cor_sig3$PSVR, plt_cor_sig3$Contribution)
ggplot(plt_cor_sig3, aes(100*Contribution, PSVR)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Signature 3 contribution (%)", y="PSVR output") +
  scale_y_continuous(limits=c(0,1.1)) +
  annotate("text", 52, 0.18, label="PCC (95% CI) = 0.10 (-0.09, 0.29) \n P = 0.3", size=5) +
  myScatterTheme

## Mutation rates per megabase:
plt_mutrates <- merge(FINAL_PREDS_TCGA, KANDOTH_MUTS, by=KEY)
t.test(Mutation.Rate...Mbp. ~ PSVR_strata, data=plt_mutrates)
ggplot(plt_mutrates, aes(PSVR_strata, Mutation.Rate...Mbp.)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA) +
  geom_jitter(width=0.25) +
  annotate("text", 1.5, 5, label="** P = 0.003", size=8) +
  labs(x="Predicted HRD Stratum (PSVR)", y="Mutation rate per megabase") +
  myBoxplotTheme

cor.test(plt_mutrates$PSVR, plt_mutrates$Mutation.Rate...Mbp.)
ggplot(plt_mutrates, aes(Mutation.Rate...Mbp., PSVR)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="Mutation rate per megabase", y="PSVR output") +
  annotate("text", 4, 0.3, label="PCC (95% CI) = 0.43 (0.22, 0.61) \n P = 2.2e-3 ***", size=5) +
  scale_y_continuous(limits=c(0,1.2), breaks=seq(0,0.9, 0.3)) +
  myScatterTheme

# ---------------------------------- Overall & Progression-free Survival  ---------------------------------- 
plt_clin_surv <- merge(FINAL_PREDS_TCGA[,c(KEY,"PSVR","PSVR_strata")], TCGA_SURV, by=KEY)

plt_clin_surv$OS[plt_clin_surv$OS.time >= ADMIN_CENSOR_Y] <- FALSE
plt_clin_surv$OS.time[plt_clin_surv$OS.time >= ADMIN_CENSOR_Y] <- ADMIN_CENSOR_Y

plt_clin_surv$PFI[plt_clin_surv$PFI.time >= ADMIN_CENSOR_Y] <- FALSE
plt_clin_surv$PFI.time[plt_clin_surv$PFI.time >= ADMIN_CENSOR_Y] <- ADMIN_CENSOR_Y

## Kaplan-Meier analysis / visualization:
kmMods <- list()
kmMods[["OS"]] <- survfit(Surv(time=OS.time, event=OS) ~ PSVR_strata, data=plt_clin_surv)
kmMods[["PFS"]] <- survfit(Surv(time=PFI.time, event=PFI) ~ PSVR_strata, data=plt_clin_surv)

for(mod in kmMods) {
  print(ggsurvplot(
    mod,
    data = plt_clin_surv,
    tables.theme = theme_cleantable(),
    
    ## Plot:
    pval = TRUE,
    pval.size = 8,
    conf.int = FALSE,
    conf.int.style = "step",
    censor = TRUE,
    size = 1.5,
    alpha = 0.7,
    
    ## Labels:
    legend = c(0.7, 0.2), #relative position
    legend.title = "Stratum",
    
    ## Global elements:
    ggtheme = mySurvTheme,
    palette = "Set1"
  ))
}
