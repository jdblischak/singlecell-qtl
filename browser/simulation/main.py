import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import os.path
import numpy as np
import pandas as pd
import scipy.stats as st
import scipy.special as sp
import sqlite3

db = '/project2/mstephens/aksarkar/projects/singlecell-qtl/browser/browser.db'

sim_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=[
  'log_mu',
  'log_mu_hat',
  'log_phi',
  'log_phi_hat',
  'logodds',
  'logodds_hat',
  'mean',
  'true_mean',
  'num_mols',
  'num_samples',
  'trial',
  'var',
  'true_var',
  'fano',
  'true_fano',
  'mean_log_cpm',
  'var_log_cpm',
]))

theoretical_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=[
  'log_mu0',
  'log_mu1',
  'log_phi0',
  'log_phi1',
  'logodds0',
  'logodds1',
  'true_mean0',
  'true_mean1',
  'true_var0',
  'true_var1',
  'true_fano0',
  'true_fano1',
]))

# This needs to be separate because it has different dimension
cpm_vs_mu_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=[
  'log_mu',
  'mean',
  'var',
]))
cpm_vs_phi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=[
  'log_phi',
  'var',
]))

num_samples_slider = bokeh.models.widgets.Slider(title='Number of samples', value=100, start=100, end=1000, step=225)
num_mols_slider = bokeh.models.widgets.Slider(title='Number of molecules', value=100e3, start=100e3, end=1000e3, step=225e3)
log_mu_slider = bokeh.models.widgets.RangeSlider(title='log(μ)', value=(-12, -6), start=-12, end=-6, step=1)
log_phi_slider = bokeh.models.widgets.RangeSlider(title='log(φ)', value=(-6, -6), start=-6, end=0, step=1)
logodds_slider = bokeh.models.widgets.RangeSlider(title='logit(π)', value=(-3, -3), start=-3, end=3, step=1)
controls = [num_samples_slider, num_mols_slider, log_mu_slider, log_phi_slider, logodds_slider]

def update(attr, old, new):
  global sim_data
  args = [x.value for x in controls[:2]] + [round(endpoint) for x in controls[2:] for endpoint in x.value]
  with sqlite3.connect(db) as conn:
    params = pd.read_sql(
      """select * from simulation 
      where num_samples == ? and num_mols == ? and 
      log_mu >= ? and log_mu <= ? and 
      log_phi >= ? and log_phi <= ? and 
      logodds >= ? and logodds <= ?""",
      conn,
      params=args)
    params['true_mean'] = params['num_mols'] * np.exp(params['log_mu'])
    params['true_var'] = params['true_mean'] + np.square(params['true_mean']) * np.exp(params['log_phi'])
    params['fano'] = params['var'] / params['mean']
    params['true_fano'] = params['true_var'] / params['true_mean']
    sim_data.data = bokeh.models.ColumnDataSource.from_df(params)
    theoretical_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(
      {
        'log_mu0': params['log_mu'].min(),
        'log_mu1': params['log_mu'].max(),
        'log_phi0': params['log_phi'].min(),
        'log_phi1': params['log_phi'].max(),
        'logodds0': params['logodds'].min(),
        'logodds1': params['logodds'].max(),
        'true_mean0': params['true_mean'].min(),
        'true_mean1': params['true_mean'].max(),
        'true_var0': params['true_var'].min(),
        'true_var1': params['true_var'].max(),
        'true_fano0': params['true_fano'].min(),
        'true_fano1': params['true_fano'].max(),
      }, index=[0]))

    eps = .5 / params['num_mols'].iloc[0]
    mean_grid = np.linspace(-12, -6, 100)
    mu = params['num_mols'].iloc[0] * np.exp(mean_grid)
    cpm_vs_mu_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame({
      'log_mu': mean_grid,
      'mean': np.log(mu + eps) - np.log(params['num_mols'].iloc[0] + 2 * eps) + 6 * np.log(10),
      'var': 
    }))

for c in controls:
  c.on_change('value', update)

# Initialize the view here
update(None, None, None)

tools = []

log_mu = bokeh.plotting.figure(width=300, height=300, tools=tools)
log_mu.scatter(source=sim_data, x='log_mu', y='log_mu_hat', color='black', size=6)
log_mu.segment(source=theoretical_data, x0='log_mu0', y0='log_mu0', x1='log_mu1', y1='log_mu1', line_width=1)
log_mu.xaxis.axis_label = 'True log(μ)'
log_mu.yaxis.axis_label = 'Estimated log(μ)'

log_phi = bokeh.plotting.figure(width=300, height=300, tools=tools)
log_phi.scatter(source=sim_data, x='log_phi', y='log_phi_hat', color='black', size=6)
log_phi.segment(source=theoretical_data, x0='log_phi0', y0='log_phi0', x1='log_phi1', y1='log_phi1', line_width=1)
log_phi.xaxis.axis_label = 'True log(φ)'
log_phi.yaxis.axis_label = 'Estimated log(φ)'

logodds = bokeh.plotting.figure(width=300, height=300, tools=tools)
logodds.scatter(source=sim_data, x='logodds', y='logodds_hat', color='black', size=6)
logodds.segment(source=theoretical_data, x0='logodds0', y0='logodds0', x1='logodds1', y1='logodds1', line_width=1)
logodds.xaxis.axis_label = 'True logit(π)'
logodds.yaxis.axis_label = 'Estimated logit(π)'

mean = bokeh.plotting.figure(width=300, height=300, tools=tools)
mean.scatter(source=sim_data, x='true_mean', y='mean', color='black', size=6)
mean.segment(source=theoretical_data, x0='true_mean0', y0='true_mean0', x1='true_mean1', y1='true_mean1', line_width=1)
mean.xaxis.axis_label = 'True mean'
mean.yaxis.axis_label = 'Estimated mean'

var = bokeh.plotting.figure(width=300, height=300, tools=tools)
var.scatter(source=sim_data, x='true_var', y='var', color='black', size=6)
var.segment(source=theoretical_data, x0='true_var0', y0='true_var0', x1='true_var1', y1='true_var1', line_width=1)
var.xaxis.axis_label = 'True variance'
var.yaxis.axis_label = 'Estimated variance'

fano = bokeh.plotting.figure(width=300, height=300, tools=tools)
fano.scatter(source=sim_data, x='true_fano', y='fano', color='black', size=6)
fano.segment(source=theoretical_data, x0='true_fano0', y0='true_fano0', x1='true_fano1', y1='true_fano1', line_width=1)
fano.xaxis.axis_label = 'True Fano factor'
fano.yaxis.axis_label = 'Estimated Fano factor'

mean_log_cpm = bokeh.plotting.figure(width=300, height=300, tools=tools)
mean_log_cpm.scatter(source=sim_data, x='log_mu', y='mean_log_cpm', color='black', size=6)
mean_log_cpm.line(source=cpm_vs_mu_data, x='log_mu', y='mean')
mean_log_cpm.xaxis.axis_label = 'True log(μ)'
mean_log_cpm.yaxis.axis_label = 'Mean log CPM'

var_log_cpm = bokeh.plotting.figure(width=300, height=300, tools=tools)
var_log_cpm.scatter(source=sim_data, x='log_mu', y='var_log_cpm', color='black', size=6)
var_log_cpm.line(source=cpm_vs_mu_data, x='log_mu', y='var')
var_log_cpm.xaxis.axis_label = 'True log(μ)'
var_log_cpm.yaxis.axis_label = 'Variance log CPM'

var_log_cpm2 = bokeh.plotting.figure(width=300, height=300, tools=tools)
var_log_cpm2.scatter(source=sim_data, x='log_phi', y='var_log_cpm', color='black', size=6)
var_log_cpm2.xaxis.axis_label = 'True log(φ)'
var_log_cpm2.yaxis.axis_label = 'Variance log CPM'

zinb_params = bokeh.models.widgets.Panel(
  child=bokeh.layouts.gridplot([[log_mu, log_phi, logodds]]),
  title='ZINB parameters')
zinb_phenos = bokeh.models.widgets.Panel(
  child=bokeh.layouts.gridplot([[mean, var, fano]]),
  title='ZINB phenotypes')
log_cpm = bokeh.models.widgets.Panel(
  child=bokeh.layouts.gridplot([[mean_log_cpm, var_log_cpm, var_log_cpm2]]),
  title='log CPM')

layout = bokeh.layouts.layout([[
  bokeh.layouts.widgetbox(*controls, width=300),
  bokeh.layouts.widgetbox(
    bokeh.models.widgets.Tabs(tabs=[zinb_params, zinb_phenos, log_cpm]),
    width=900),
]])

doc = bokeh.io.curdoc()
doc.title = 'Simulation results'
doc.add_root(layout)
