# Fraction of Genome Altered (FGA)
# Copyright (C) 2019-2024 Y. David Chen & Christensen Lab. All rights reserved.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(reshape2)
source("../utils.R")

source("../plotThemes.R")
myScatterTheme$panel.spacing.x <- unit(2, "lines")
myBoxplotTheme$axis.title.x <- element_text(size=20,color="black")

run_corr_fit <- function(x, y) {
  fit_lm <- lm(x ~ y)
  print(summary(fit_lm))
  print("---------------------------------")
  print(cor.test(x, y, method="pearson"))
}

KEY <- "Sample_ID"
COLS_FGA <- c("nseg_raw","fga_raw","nseg_auto","fga_auto")
COLS_PLT <- c("Sample_ID","fga_auto","HRD","MLPA.mean","RPMM")

# ------------------ Load FGA & Clinical Data ------------------
## Load Cohort Annotations:
load_cohort1(); load_cohort2()
rm(cohort1_betas, cohort2_betas)

## Update cohort annotations with FGA:
fga_cohort1 <- read.csv(paste0(DIR_CNA,"cohort1_fga_stats.csv"), row.names=1, stringsAsFactors=FALSE)
fga_cohort1[KEY] <- rownames(fga_cohort1)
cohort1_covars <- merge(cohort1_covars, fga_cohort1[,c(KEY,COLS_FGA)], by=KEY)

fga_cohort2 <- read.csv(paste0(DIR_CNA,"cohort2_fga_stats.csv"), row.names=1, stringsAsFactors=FALSE)
fga_cohort2[KEY] <- rownames(fga_cohort2)
cohort2_covars[KEY] <- cohort2_covars$Tube.number
cohort2_covars <- merge(cohort2_covars, fga_cohort2[,c(KEY,COLS_FGA)], by=KEY)

# ------------------ Data visualization ------------------
cohorts_fga <- rbind(
  cbind(cohort1_covars[ , COLS_PLT], Cohort="Cohort 1"),
  cbind(cohort2_covars[ , COLS_PLT], Cohort="Cohort 2")
)

## By HRD: 
ggplot(cohorts_fga, aes(x=Cohort, y=fga_auto, color=HRD)) +
  geom_boxplot(outlier.colour=NA, outlier.size=NA) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.25)) +
  geom_smooth(method="lm") +
  scale_color_brewer(palette="Set1") +
  ylab("Fraction of Genome Altered") +
  myBoxplotTheme

## By RPMM Clusters:
ggplot(cohorts_fga, aes(RPMM, fga_auto)) +
  geom_boxplot(outlier.colour=NA, outlier.size=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort, scales="free_x") +
  ylab("Fraction of Genome Altered") +
  myBoxplotTheme

#------------------ Statistical Tests ------------------
t.test(fga_auto ~ HRD, cohort1_covars)
t.test(fga_auto ~ HRD, cohort2_covars)

run_corr_fit(cohort2_covars$MLPA.mean, cohort2_covars$fga_auto)
run_corr_fit(cohort1_covars$MLPA.mean, cohort1_covars$fga_auto)

summary( aov(fga_auto ~ RPMM, data=cohort1_covars) )
summary( aov(fga_auto ~ RPMM, data=cohort2_covars) )

# ------------------ Composite Figures for Manuscript ------------------
plt_cohorts <- rbind(
  cbind(cohort1_covars[ , c("fga_auto","HRD","MLPA.mean")], Cohort="Cohort 1"),
  cbind(cohort2_covars[ , c("fga_auto","HRD","MLPA.mean")], Cohort="Cohort 2")
)

# png("~/Downloads/fga_vs_hrd.png", width=15, height=7.5, units="in", res=300);
ggplot(plt_cohorts, aes(MLPA.mean, fga_auto)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="HRD Probability", y="Fraction of Genome Altered") +
  geom_text(data=data.frame(x=0.3,y=0.6,Cohort="Cohort 1",lab="PCC (95% CI) = 0.70 (0.46, 0.84) \n P = 9.50e-6"), aes(x=x,y=y,label=lab), size=5) +
  geom_text(data=data.frame(x=0.3,y=0.6,Cohort="Cohort 2",lab="PCC (95% CI) = 0.50 (0.28, 0.67) \n P = 5.58e-5"), aes(x=x,y=y,label=lab), size=5) +
  facet_wrap( ~ Cohort, ncol=2) +
  myScatterTheme
# dev.off()

# png("~/Downloads/fga_by_hrd_strata.png", width=5, height=10, units="in", res=300)
ggplot(plt_cohorts, aes(HRD, fga_auto)) +
  geom_boxplot(outlier.colour=NA, outlier.size=NA) +
  geom_jitter(width=0.25) +
  geom_text(data=data.frame(x=1.5,y=0.7,Cohort="Cohort 1",lab="***"), aes(x=x,y=y,label=lab), size=12, fontface=2) +
  geom_text(data=data.frame(x=1.5,y=0.7,Cohort="Cohort 2",lab="**"), aes(x=x,y=y,label=lab), size=12, fontface=2) +
  ylab("Fraction of Genome Altered") +
  scale_color_brewer(palette="Set1") +
  facet_wrap(~ Cohort) +
  myBoxplotTheme
# dev.off()
