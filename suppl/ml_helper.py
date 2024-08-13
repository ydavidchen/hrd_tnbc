__title__ = "Python Machine Learning Helper Functions for HRD"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import logit, expit

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR

## Helper functions:
def z_norm_by_col(X):
    """
    Standard-normalize by column mean & s.d.
    """
    X_mean = X.mean(axis=0)
    X_stdev = X.std(axis=0)
    return (X - X_mean) / X_stdev

def draw_scatter(ypred, ytrue, sz=7, title=None):
    plt.figure(figsize=[sz,sz], dpi=200)
    plt.scatter(ypred, ytrue, alpha=0.6, color="black")
    if title is not None:
        plt.title(title, size=15)
    plt.xlabel("Prediction", size=15)
    plt.ylabel("Actual", size=15)
    plt.show()

def show_grid_stats(clf):
    print("Average score on validatoin folds for each hyperparameter set:")
    means, stds, params = clf.cv_results_['mean_test_score'], clf.cv_results_['std_test_score'], clf.cv_results_['params']
    for mean, std, param in zip(means, stds, params):
        print("%0.3f (+/-%0.03f) for %r" % (mean, std * 2, param))
    print()
    print("The best set of hyperparameters are:" + str(clf.best_params_))
    print("...with the following scores: %0.4f" % clf.best_score_)

def evaluate_regressor(preds, truths, plotTitle=None, fontsize=13, out_prefix=None):
    from scipy.stats import pearsonr;
    from sklearn.metrics import r2_score, explained_variance_score, mean_squared_error;

    def apply_constraint(val, upper=1, lower=0):
        if val > upper:
            return upper;
        if val < lower:
            return lower;
        return val;

    ## Compute & display correlation statistics
    rmse = np.sqrt(mean_squared_error(truths, preds))
    pcc, _ = pearsonr(truths, preds)
    m, b = np.polyfit(preds, truths, 1)
    rsq = r2_score(truths, preds)
    var_exp = explained_variance_score(truths, preds)

    if plotTitle is not None:
        plt.figure(figsize=(6,6), dpi=100)
        plt.scatter(preds, truths, c='b', alpha=0.5)
        plt.xlabel("Predicted Value", fontsize=fontsize)
        plt.ylabel("Actual Value", fontsize=fontsize)

        ## Trend lines:
        axes = plt.gca()
        X_plot = np.linspace(axes.get_xlim()[0], axes.get_xlim()[1], 100)
        plt.plot(X_plot, m * X_plot + b, 'b--', label='Model')
        plt.plot(X_plot, X_plot, 'k--', alpha=0.5, label='Perfect fit')
        plt.legend(loc='lower right')

        ## Text labels:
        metric_text = "\n RMSE = %.3f, Pearson Cor = %.3f, R-squared = %.3f" % (rmse, pcc, rsq)
        plt.title(plotTitle+metric_text, fontsize=fontsize)

        ## Optional: Save figure to file
        if out_prefix is not None:
            plt.savefig(out_prefix + plotTitle + ".png")

        plt.show()
        plt.close()

    return {'PCC':pcc, 
            'RMSE':rmse, 
            'R2':apply_constraint(rsq),
            'varExpl':apply_constraint(var_exp), 
            'slope':m, 
            'intercept':b}

def parse_by_dataset(full_data, sub_annot, label_col=None, transform_label=False):
    """ Wrapper function to assign methylation data & associated label """
    sub_data = full_data[full_data.index.isin(sub_annot.index)];
    if label_col is None:
        return sub_data;

    if all(sub_annot.index == sub_data.index):
        print("Condition met! Proceed to extracting regression labels...")
        y = sub_annot[label_col]
        if transform_label:
            print("Logit-transforming labels...")
            y = logit(y)
        y.index = sub_annot.index;
    else:
        raise ValueError("Condition failed! Requires matching and/or subsetting first")
    return sub_data, y;


def build_and_test_regressor(regressor, Xtrain, Xtest, ytrain, ytest, plotTitle=None):
    """ Applies an initialized sklearn regressor to training & test data """
    regressor.fit(Xtrain, ytrain)
    test_pred = expit(regressor.predict(Xtest).ravel())
    test_groundtruth = expit(ytest.ravel())
    return evaluate_regressor(test_pred, test_groundtruth, plotTitle)

def ab_test_wrapper(regressorA, regressorB, X, y, test_size=30, n_iter=30):
    """ Performs A/B test with same train-test split """
    res_A, res_B = [], []
    for n in range(n_iter):
        Xtrain_k, Xtest_k, ytrain_k, ytest_k = train_test_split(X, y, test_size=test_size, random_state=None)
        res_A.append(build_and_test_regressor(regressorA, Xtrain_k, Xtest_k, ytrain_k, ytest_k))
        res_B.append(build_and_test_regressor(regressorB, Xtrain_k, Xtest_k, ytrain_k, ytest_k))

    return pd.DataFrame(res_A), pd.DataFrame(res_B)
