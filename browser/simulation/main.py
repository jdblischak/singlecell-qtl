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
  "log_mu",
  "log_mu_hat",
  "log_phi",
  "log_phi_hat",
  "logodds",
  "logodds_hat",
  "mean",
  "num_mols",
  "num_samples",
  "trial",
  "var",
]))

num_samples_slider = bokeh.models.widgets.Slider(title="Number of samples", value=100, start=100, end=1000, step=225)
num_mols_slider = bokeh.models.widgets.Slider(title="Number of molecules", value=100e3, start=100e3, end=1000e3, step=225e3)
log_mu_slider = bokeh.models.widgets.RangeSlider(title="log(μ)", value=(-12, -6), start=-12, end=-6, step=1)
log_phi_slider = bokeh.models.widgets.RangeSlider(title="log(φ)", value=(-6, -6), start=-6, end=0, step=1)
logodds_slider = bokeh.models.widgets.RangeSlider(title="logit(π)", value=(-3, -3), start=-3, end=3, step=1)
controls = [num_samples_slider, num_mols_slider, log_mu_slider, log_phi_slider, logodds_slider]

def update(attr, old, new):
  global sim_data
  args = [x.value for x in controls[:2]] + [endpoint for x in controls[2:] for endpoint in x.value]
  with sqlite3.connect(db) as conn:
    sim_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      """select * from simulation 
      where num_samples == ? and num_mols == ? and 
      log_mu >= ? and log_mu <= ? and 
      log_phi >= ? and log_phi <= ? and 
      logodds >= ? and logodds <= ?""",
      conn,
      params=args))

for c in controls:
  c.on_change('value', update)

# Initialize the view here
update(None, None, None)

tools = []

log_mu = bokeh.plotting.figure(width=300, height=300, tools=tools)
log_mu.scatter(source=sim_data, x='log_mu', y='log_mu_hat', color='black', size=6)
log_mu.segment(x0=-12, y0=-12, x1=-6, y1=-6, color='red', line_width=1)
log_mu.xaxis.axis_label = 'True log(μ)'
log_mu.yaxis.axis_label = 'Estimated log(μ)'

log_phi = bokeh.plotting.figure(width=300, height=300, tools=tools)
log_phi.scatter(source=sim_data, x='log_phi', y='log_phi_hat', color='black', size=6)
log_phi.segment(x0=-6, y0=-6, x1=0, y1=0, color='red', line_width=1)
log_phi.xaxis.axis_label = 'True log(φ)'
log_phi.yaxis.axis_label = 'Estimated log(φ)'

logodds = bokeh.plotting.figure(width=300, height=300, tools=tools)
logodds.scatter(source=sim_data, x='logodds', y='logodds_hat', color='black', size=6)
logodds.segment(x0=-3, y0=-3, x1=3, y1=3, color='red', line_width=1)
logodds.xaxis.axis_label = 'True logit(π)'
logodds.yaxis.axis_label = 'Estimated logit(π)'

layout = bokeh.layouts.layout([[
  bokeh.layouts.widgetbox(*controls, width=300),
  bokeh.layouts.gridplot([[log_mu, log_phi, logodds]]),
]])

doc = bokeh.io.curdoc()
doc.title = 'Simulation results'
doc.add_root(layout)
