# Supervised Machine Learning of Subject Race based on 59 SNP Probes

import pandas as pd
import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import ConfusionMatrixDisplay, confusion_matrix, classification_report

def evaluate_clf(clf):
    clf.fit(Xtrain, ytrain)
    yTestPred = clf.predict(Xtest)
    print(classification_report(ytest, yTestPred))
    ConfusionMatrixDisplay(confusion_matrix(ytest, yTestPred, labels=clf.classes_)).plot()
    
def predict_unk(clf, X_unk):
    clf.fit(X, encoder.transform(y))
    ypred = encoder.inverse_transform(clf.predict(X_unk))
    return(pd.Series(ypred, index=X_unk.index))

DIR_SNP = "********* MASKED ********* "

## Load & Split Data:
X = pd.read_csv(DIR_SNP+"sklearn_tcga_snps.csv", index_col=0)
y = X.pop("Race.Category")

X.head()

encoder = LabelEncoder()
encoder.fit(y)
encoder.classes_

Xtrain, Xtest, ytrain, ytest = train_test_split(X, encoder.transform(y), test_size=0.2)


## Elastic Net
cv_lr = GridSearchCV(
    LogisticRegression(), 
    param_grid = {'penalty': ['l1','l2'], 'C': [0.001,0.01,0.1,1,10,100]}, 
    cv = 5
)
cv_lr.fit(Xtrain, ytrain)
cv_lr.best_params_, cv_lr.best_score_

evaluate_clf(cv_lr.best_estimator_)


## Random Forest
cv_rf = GridSearchCV(
    RandomForestClassifier(),
    param_grid = {
        'max_depth': np.arange(10, 100, 10),
        'n_estimators': [25, 50, 100]
    },
    cv = 5
)

cv_rf.fit(Xtrain, ytrain)
cv_rf.best_params_, cv_rf.best_score_

evaluate_clf(cv_rf.best_estimator_)


## SVM 
cv_svm = GridSearchCV(
    SVC(),
    param_grid =  [
      {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
      {'C': [1, 10, 100, 1000], 'gamma': [1, 0.01, 0.001, 0.0001], 'kernel': ['rbf']},
    ],
    cv = 5
)

cv_svm.fit(Xtrain, ytrain)
cv_svm.best_params_, cv_svm.best_score_

evaluate_clf(cv_svm.best_estimator_)


## Make Predictions on Unknown Samples
snp_cohort1 = pd.read_csv(DIR_SNP+"sklearn_cohort1_snps.csv", index_col=0)
assert all(snp_cohort1.columns == X.columns) #checkpoint
snp_cohort1.head()

snp_cohort2 = pd.read_csv(DIR_SNP+"sklearn_cohort2_snps.csv", index_col=0)
assert all(snp_cohort2.columns == X.columns) #checkpoint
snp_cohort2.head()

res_cohort1 = pd.DataFrame({
    "Elastic": predict_unk(cv_lr.best_estimator_, snp_cohort1),
    "RandomForest": predict_unk(cv_rf.best_estimator_, snp_cohort1),
    "SVM": predict_unk(cv_svm.best_estimator_, snp_cohort1)
})

res_cohort1

res_cohort2 = pd.DataFrame({
    "Elastic": predict_unk(cv_lr.best_estimator_, snp_cohort2),
    "RandomForest": predict_unk(cv_rf.best_estimator_, snp_cohort2),
    "SVM": predict_unk(cv_svm.best_estimator_, snp_cohort2)
})

res_cohort2 
