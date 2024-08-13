# In-depth Evaluation HRD-Polynomial Support Vector Regressor (PSVR) & Make Predictions on TCGA
# Objectives:
# 1. A/B testing of best model (PSVR) against baseline (OLS)
# 2. Consistency across 100 iterations of random train-test splits
# 3. Re-train final model and make regression predictions on ALL Cohorts 1-2 tumors & TCGA. Save exported scores to file

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import logit, expit
from sklearn.model_selection import train_test_split, GridSearchCV

from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR

from ml_helper import *

PATH_MVALS = "********* MASKED *********"
PATH_CLIN = "********* MASKED *********"
DIR_OUT = "**************** MASKED ****************"
SKEY = ['Mean', 'Median', 'SD']

## Initialize selected model & baseline:
SELECTED_PSVR = FINAL_PSVR = SVR(kernel='poly', degree=5, C=1, gamma=1e-4)
BASELINE = FINAL_BASE = LinearRegression()

## Load & parse out labeled data: 
full_data = pd.read_csv(PATH_MVALS, index_col=0)
full_data.index = full_data.index.map(str)

full_annot = pd.read_csv(PATH_CLIN, index_col=0)
full_annot.index = full_annot.index.map(str)

trainval_annot = full_annot[full_annot.Dataset != "TCGA-TNBCs"]
trainval_annot.head()

X, y = parse_by_dataset(full_data, trainval_annot, "MLPA", True)

## Parse out TCGA as the unknown for HRD high/low assignment:
tcga_annot = full_annot[full_annot.Dataset=="TCGA-TNBCs"]
X_unknown = parse_by_dataset(full_data, tcga_annot)

print("Training with %d samples and %d features" % Xtrain.shape)
print("Evaluating with %d samples and %d features" % Xtest.shape)
print("Making predictions on %d TCGA samples and %d features" % X_unknown.shape)

## Tasks 1 & 2: A/B Testing with 100-trial Stability Evaluation
tmp_psvr, tmp_ctrl = ab_test_wrapper(SELECTED_PSVR, LinearRegression(), X, y, 100)

pd.concat(
    [tmp_psvr.mean(axis=0), tmp_psvr.median(axis=0), tmp_psvr.std(axis=0)], 
    axis = 1, 
    keys = SKEY
)

pd.concat(
    [tmp_ctrl.mean(axis=0), tmp_ctrl.median(axis=0), tmp_ctrl.std(axis=0)], 
    axis = 1, 
    keys = SKEY
)

# tmp_psvr.to_csv(DIR_OUT+"consistency_eval_psvr.csv")
# tmp_ctrl.to_csv(DIR_OUT+"consistency_eval_baseline.csv")

## Task 3: Retrain final models w/ all labeled TNBCs and make predictions on TCGA unknowns
## Sanity check:
build_and_test_regressor(FINAL_PSVR, X, X, y, y, "PSVR")
build_and_test_regressor(BASELINE, X, X, y, y, "Baseline")

## Final modeling:
FINAL_PSVR.fit(X, y)
FINAL_BASE.fit(X, y)

## Make predictions on TCGA & export:
TCGA_PRED_PSVR = expit(FINAL_PSVR.predict(X_unknown).ravel())
TCGA_PRED_BASE = expit(FINAL_BASE.predict(X_unknown).ravel())

df_export = pd.concat(
    [pd.Series(X_unknown.index.values), pd.Series(TCGA_PRED_PSVR), pd.Series(TCGA_PRED_BASE)], 
    keys = ["bcr_patient_barcode", "PSVR", "Baseline"],
    axis = 1
)
# df_export.to_csv(DIR_OUT+"TCGA_predictions.csv") #for downstream R analysis
