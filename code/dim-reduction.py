"""Script to reproduce Fig. 1B and 1D

Usage: 

source activate scqtl
MPLBACKEND=Agg python dim-reduction.py

Output:

Files figure-1b.png, figure-1d.png

Author: Abhishek Sarkar <aksarkar@uchicago.edu>

"""

import colorcet
import functools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd
import scipy.stats as st
import sklearn.decomposition as skd
import sklearn.linear_model as sklm

def extract_covars(annotations):
  """Return a Series of covariate values

  Map keys of annotations to human-readable names for the figure

  """
  return pd.Series({
    'Reads': annotations['umi'],
    'Molecules': annotations['molecules'],
    'Mapped %': annotations['umi'] / annotations['mapped'],
    'Genes detected': annotations['detect_hs'],
  })

def correlation(pcs, cont_covars):
  """Return squared correlation between principal components and covariates

  pcs - DataFrame (n x k)
  cont_covars - DataFrame (n x q)

  """
  result = []
  for i in pcs:
    for j in cont_covars:
      keep = np.isfinite(cont_covars[j].values)
      result.append([i, j, np.square(st.pearsonr(pcs[i][keep], cont_covars[j][keep]))[0]])
  return pd.DataFrame(result, columns=['pc', 'covar', 'corr'])

def categorical_r2(loadings, annotations, key, name):
  """Return coefficient of determination of regression of PC loadings against
onehot-encoded categorical variable

  """
  categories = sorted(annotations[key].unique())
  onehot = np.zeros((annotations.shape[0], len(categories)), dtype=np.float32)
  onehot[np.arange(onehot.shape[0]), annotations[key].apply(lambda x: categories.index(x))] = 1
  m = sklm.LinearRegression(fit_intercept=True, copy_X=True).fit(onehot, loadings)
  return pd.DataFrame({
      'pc': np.arange(10),
      'covar': name,
      'corr': 1 - np.square(loadings - m.predict(onehot)).sum(axis=0) / np.square(loadings - loadings.mean(axis=0)).sum(axis=0)})

def sample_feature_means(obs, max_iters=10):
  """Fit per-feature and per-sample means

  x_ij ~ N(u_i + v_j, sigma^2)

  Inputs:

  x - ndarray (num_samples, num_features)

  Returns:
  sample_means - ndarray (num_samples, 1)
  feature_means - ndarray (num_features, 1)

  """
  n, p = obs.shape
  sample_means = np.zeros((n, 1))
  feature_means = np.zeros((1, p))
  resid = obs.var()
  llik = -.5 * np.sum(np.log(2 * np.pi * resid) + (obs - feature_means - sample_means) / resid)
  # Coordinate ascent on llik
  for _ in range(max_iters):
    sample_means = np.mean(obs - feature_means, axis=1, keepdims=True)
    feature_means = np.mean(obs - sample_means, axis=0, keepdims=True)
    # By default, np divides by N, not N - 1
    resid = np.var(obs - feature_means - sample_means)
    update = -.5 * np.sum(np.log(2 * np.pi * resid) + (obs - feature_means - sample_means) / resid)
    print(update)
    if np.isclose(update, llik):
      return sample_means, feature_means.T
    else:
      llik = update
  raise ValueError("Failed to converge")

def plot_pca_covar_corr(pca, corr):
  """Plot heatmap of correlation between PCs and covariates and PVE of PCs"""
  plt.clf()
  fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [.25, .75]}, sharex=True)
  fig.set_size_inches(4, 5)
  ax[0].bar(np.arange(len(pca.components_)), pca.explained_variance_ratio_)
  ax[0].set_xticks(np.arange(len(pca.components_)))
  ax[0].set_xticklabels([str(x) for x in np.arange(1, len(pca.components_) + 1)])
  ax[0].set_xlabel('Principal component')
  ax[0].set_ylabel('PVE')

  im = ax[1].imshow(corr.values, cmap=colorcet.cm['fire'], vmin=0, vmax=1, aspect='auto')
  cb = plt.colorbar(im, ax=ax[1], orientation='horizontal')
  cb.set_label('Squared correlation')
  ax[1].set_xlabel('Principal component')
  ax[1].set_yticks(np.arange(corr.shape[0]))
  ax[1].set_yticklabels(corr.index)
  ax[1].set_ylabel('Covariate')

  plt.gcf().tight_layout()

def plot_bicentered_pca(log_cpm, annotations, cont_covars, cat_covars):
  """Bicenter the data, then plot the PCA correlation + PVE"""
  sample_means, feature_means = sample_feature_means(log_cpm.values.T)
  ppca = skd.PCA(n_components=10)
  loadings = ppca.fit_transform(log_cpm.values.T - sample_means - feature_means.T)
  corr = pd.concat(
    [correlation(pd.DataFrame(loadings), cont_covars)] +
    [categorical_r2(loadings, annotations, k, n) for k, n in cat_covars])
  corr = corr.pivot(index='covar', columns='pc')
  plot_pca_covar_corr(ppca, corr)

if __name__ == '__main__':
  # Get the root of the singlecell-qtl repository
  basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')

  # Read the data
  umi = pd.read_table('{}/data/scqtl-counts.txt.gz'.format(basedir), index_col=0)
  annotations = pd.read_table('{}/data/scqtl-annotation.txt'.format(basedir))
  keep_samples = pd.read_table('{}/data/quality-single-cells.txt'.format(basedir), index_col=0, header=None)
  keep_genes = pd.read_table('{}/data/genes-pass-filter.txt'.format(basedir), index_col=0, header=None)
  umi = umi.loc[keep_genes.values.ravel(),keep_samples.values.ravel()]
  annotations = annotations.loc[keep_samples.values.ravel()]

  # Normalize the data
  libsize = annotations['mol_hs'].values
  pseudocount = .5 * libsize / libsize.mean()
  log_cpm = (np.log(umi + pseudocount) - np.log(libsize + 2 * pseudocount) + 6 * np.log(10)) / np.log(2)

  # Extract the covariates of intereset
  cont_covars = annotations.apply(extract_covars, axis=1)
  cat_covars = list(zip(annotations[['batch', 'experiment', 'chip_id', 'well']],
                        ['Batch', 'C1 chip', 'Individual', 'Well']))

  # Reproduce Fig. 1B
  ppca = skd.PCA(n_components=10)
  loadings = ppca.fit_transform(log_cpm.values.T)
  corr = pd.concat(
    [correlation(pd.DataFrame(loadings), cont_covars)] +
    [categorical_r2(loadings, annotations, k, n) for k, n in cat_covars])
  corr = corr.pivot(index='covar', columns='pc')
  plot_pca_covar_corr(ppca, corr)
  plt.savefig('figure-1b.png', dpi=144)

  # Reproduce Fig. 1D
  normalized_percentiles = np.percentile(log_cpm, 100 * np.linspace(0, 1, 5), axis=0)
  log_cpm_residuals = log_cpm.apply(lambda y: y - sklm.LinearRegression(fit_intercept=True).fit(normalized_percentiles.T, y).predict(normalized_percentiles.T), axis=1)
  plot_bicentered_pca(log_cpm_residuals, annotations, cont_covars, cat_covars)
  plt.savefig('figure-1d.png', dpi=144)
