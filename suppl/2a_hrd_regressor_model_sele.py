# Model Selection for Methylation-based Continuous HRD Predictor

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import logit, expit

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LinearRegression, ElasticNet
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor

from ml_helper import *

## Constants:
SEED = 200117
PATH_MVALS = "********* MASKED *********"
PATH_CLIN = "********* MASKED *********"
CV_SCORE = "neg_mean_squared_error"

## Grids to search:
ENET_PARAM_GRID = {'alpha':np.arange(0.05,1.05,0.05)}
RF_PARAM_GRID = {'max_depth': np.arange(10, 100, 10), 'n_estimators': [25, 50, 100, 200]}
SVM_PARAM_GRID = [
    {'C': [1,10,100,1000], 'gamma':[0.001, 0.0001], 'kernel':['linear']},
    {'C': [1,10,100,1000], 'gamma':[0.001, 0.0001], 'degree':[2,3,4,5], 'kernel':['poly']}
]
KNN_PARAM_GRID = {'n_neighbors':range(5,20), 'metric':['euclidean','manhattan']}


## Load TNBC cohort data prepared for Python ML
full_data = pd.read_csv(PATH_MVALS, index_col=0)
full_data.index = full_data.index.map(str)

full_annot = pd.read_csv(PATH_CLIN, index_col=0)
full_annot.index = full_annot.index.map(str)

trainval_annot = full_annot[full_annot.Dataset != "TCGA-TNBCs"]
trainval_annot.head()

X, y = parse_by_dataset(full_data, trainval_annot,"MLPA", True)

if not all(y.index == X.index):
    raise ValueError("Condition failed! Requires matching and/or subsetting first")
    
Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, test_size=30, random_state=SEED)
print("Training with %d samples and %d features" % Xtrain.shape)
print("Evaluating with %d samples and %d features" % Xtest.shape)


## Unsupervised Exploration with PCA and t-SNE
HRD_cohorts = pd.Series(["Yes" if y_raw > 0.50 else "No" for y_raw in expit(y)])

X_pca = PCA(n_components=2).fit_transform(expit(X))
X_tsne = TSNE(n_components=2, perplexity=60).fit_transform(expit(X))

for X_embedded in X_pca, X_tsne:
    plt.figure(figsize=(6,6), dpi=100)
    plt.scatter(X_embedded[HRD_cohorts=="Yes",0], X_embedded[HRD_cohorts=="Yes",1], label="Yes", color='r')
    plt.scatter(X_embedded[HRD_cohorts=="No",0], X_embedded[HRD_cohorts=="No",1], label="No", color='b')
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")
    plt.legend(loc="lower right", title="HRD")
    plt.show()

## Linear Regression - baseline/neg control no CV
ols = LinearRegression()
ols.fit(Xtrain, ytrain)
evaluate_regressor(
    expit(ols.predict(Xtest).ravel()), 
    expit(ytest.ravel()), 
    "Linear Regression"
)

## Elastic-net Regression
reg_enet_grid = GridSearchCV(ElasticNet(), ENET_PARAM_GRID, cv=5, scoring=CV_SCORE)
reg_enet_grid.fit(Xtrain, ytrain)
show_grid_stats(reg_enet_grid)

reg_enet = reg_enet_grid.best_estimator_
reg_enet.fit(Xtrain, ytrain)
evaluate_regressor(
    expit(reg_enet.predict(Xtest).ravel()), 
    expit(ytest.ravel()), 
    "Elastic Net"
)

## Random Forest
reg_rf_grid = GridSearchCV(RandomForestRegressor(), RF_PARAM_GRID, cv=5, scoring=CV_SCORE)
reg_rf_grid.fit(Xtrain, ytrain)
show_grid_stats(reg_rf_grid)

reg_rf = reg_rf_grid.best_estimator_
reg_rf.fit(Xtrain, ytrain)
evaluate_regressor(
    expit(reg_rf.predict(Xtest).ravel()), 
    expit(ytest.ravel()), 
    "Random Forest"
)

## Support Vector Regression
reg_svm_grid = GridSearchCV(SVR(), SVM_PARAM_GRID, cv=5, scoring=CV_SCORE)
reg_svm_grid.fit(Xtrain, ytrain.ravel())
show_grid_stats(reg_svm_grid)

reg_svm = reg_svm_grid.best_estimator_
reg_svm.fit(Xtrain, ytrain.ravel())
evaluate_regressor(
    expit(reg_svm.predict(Xtest).ravel()), 
    expit(ytest.ravel()), 
    "Support Vector Regression"
)

## KNN Regression
reg_knn_grid = GridSearchCV(KNeighborsRegressor(), KNN_PARAM_GRID, cv=5, scoring=CV_SCORE)
reg_knn_grid.fit(Xtrain, ytrain.ravel())
show_grid_stats(reg_knn_grid)

reg_knn = reg_knn_grid.best_estimator_
reg_knn.fit(Xtrain, ytrain.ravel())
evaluate_regressor(
    expit(reg_knn.predict(Xtest).ravel()), 
    expit(ytest.ravel()), 
    "KNN"
)
