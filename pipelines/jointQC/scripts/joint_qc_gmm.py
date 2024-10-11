import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import silhouette_score
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy import linalg
#import itertools

from sklearn.mixture import BayesianGaussianMixture
from sklearn.datasets import make_classification
import numpy as np
from scipy.special import logsumexp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parse input arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = " ")
parser.add_argument('--metrics', type = str, help="")
parser.add_argument('--model', type = str, help="")
parser.add_argument('--prefilter_barcodes', type = str, help="")
parser.add_argument('--outdir', type = str, help="")

args = parser.parse_args()

METRICS_FILE = args.metrics
PREFILTER_BARCODES = args.prefilter_barcodes
OUTDIR = args.outdir
MODEL = args.model


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Gaussian Mixture models to select high quality nuclei based on qc metrics
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### GMM Functions
def gmm_bic_score(estimator, X):
    """Callable to pass to GridSearchCV that will use the BIC score."""
    # Make it negative since GridSearchCV expects a score to maximize
    return -estimator.bic(X)

def bgmm_bic_score(estimator, X):
    N, n_features = X.shape
    log_likelihood = estimator.score(X) * N  # score(X) gives avg log-likelihood per sample
    n_components = estimator.n_components
    covariance_type = estimator.covariance_type
    
    if covariance_type == 'full':
        n_parameters = n_components * n_features * (n_features + 1) / 2
    elif covariance_type == 'tied':
        n_parameters = n_features * (n_features + 1) / 2
    elif covariance_type == 'diag':
        n_parameters = n_components * n_features
    elif covariance_type == 'spherical':
        n_parameters = n_components
        
    n_parameters += n_features * n_components + n_components - 1
    
    bic = -2 * log_likelihood + n_parameters * np.log(N)
    return -bic

def plot_scores(grid_search, plotfile):
    # Collect BIC values across models
    df = pd.DataFrame(grid_search.cv_results_)[
        ["param_n_components", "param_covariance_type", "mean_test_score"]
    ]
    df["mean_test_score"] = -df["mean_test_score"]
    df = df.rename(
        columns={
            "param_n_components": "Number of components",
            "param_covariance_type": "Type of covariance",
            "mean_test_score": "BIC score",
        }
    )
    df.sort_values(by="BIC score").head()
    # Plot the BIC values
    sns.catplot(
        data=df,
        kind="bar",
        x="Number of components",
        y="BIC score",
        hue="Type of covariance",
    )
    plt.savefig(plotfile)


def plot_scores_bgmm(grid_search, plotfile):
    # Collect BIC values across models
    df = pd.DataFrame(grid_search.cv_results_)[
        ["param_n_components", "param_covariance_type", "param_weight_concentration_prior", "mean_test_score"]
    ]
    df["mean_test_score"] = -df["mean_test_score"]
    df = df.rename(
        columns={
            "param_n_components": "Number of components",
            "param_covariance_type": "Type of covariance",
            "param_weight_concentration_prior": "Weight concentration prior",
            "mean_test_score": "BIC score",
        }
    )
    df.sort_values(by="BIC score").head()
    # Plot the BIC values
    sns.catplot(
        data=df,
        kind="bar",
        x="Number of components",
        y="BIC score",
        hue="Weight concentration prior"
    )
    plt.savefig(plotfile)


def run_all(metrics, estimator, d):

    nickname = estimator + "__" + "__".join(metrics)

    # Set up input np array
    import numpy as np
    X_full = d[metrics].to_numpy()
    X_nomiss = X_full[~np.isnan(X_full).any(axis=1)]
    d_filt = d[~np.isnan(X_full).any(axis=1)]
    X = np.log10(X_nomiss+0.001)

    # Run gridsearch
    if estimator == "gmm":
        param_grid = {
            "n_components": range(1, 3),
            "covariance_type": ["spherical", "tied", "diag", "full"],
        }
        gmm = GaussianMixture()
        grid_search = GridSearchCV(
            gmm, param_grid=param_grid, scoring=gmm_bic_score
        )
        grid_search.fit(X)
        # Plot BIC scores
        plot_scores(grid_search, OUTDIR + nickname + "_bic_vals.png")
    elif estimator == "bgmm":
        param_grid = {
            "n_components": range(1, 3),
            "covariance_type": ["spherical", "tied", "diag", "full"],
            #"weight_concentration_prior_type": ['dirichlet_process'],
            "weight_concentration_prior": [0.01, 0.1, 0.5]
        }        
        bgmm = BayesianGaussianMixture(random_state=138)
        grid_search = GridSearchCV(
            bgmm, param_grid=param_grid, scoring=bgmm_bic_score
        )
        grid_search.fit(X)
        # Plot BIC scores
        plot_scores_bgmm(grid_search, OUTDIR + nickname + "_bic_vals.png")
    else:
        raise("This estimator is not supported.")
    
    # Save model predictions
    Y = grid_search.predict(X)
    d_filt.loc[:,'bestfit'] = Y.tolist()
    d_filt.to_csv(OUTDIR + nickname + "_bestfit.csv", index=False)
    


# ## Testing
# MODEL = "gmm__rna_umis__atac_hqaa"
# METRICS_FILE = "results/joint_qc_gmm/10965-VD-2/joint_qc.csv"
# PREFILTER_BARCODES = "results/droplet_utils/10965-VD-2/barcodes_nuclei.txt"
# OUTDIR = "results/joint_qc_gmm/10965-VD-2/"

# Read in data to use for modeling
d = pd.read_csv(METRICS_FILE, index_col=None)

# Filter for EmptyDrops barcodes
with open(PREFILTER_BARCODES, 'r') as file:
    barcodes_to_keep = [line.strip() for line in file]

d_filt = d[d['barcode'].isin(barcodes_to_keep)]

# Parse model
model_args = MODEL.split("__")
estimator = model_args[0]
metrics = model_args[1:]

# Run
run_all(metrics = metrics, estimator = estimator, d = d_filt)
