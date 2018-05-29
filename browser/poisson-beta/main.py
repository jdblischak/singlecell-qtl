import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import numpy as np
import pandas as pd
import scipy.stats as st

rate0_slider = bokeh.models.widgets.Slider(
  title='log10 on-off transition rate',
  value=0,
  start=-3,
  end=3,
  step=.25)
rate1_slider = bokeh.models.widgets.Slider(
  title='log10 off-on transition rate',
  value=0,
  start=-3,
  end=3,
  step=.25)
rate_slider = bokeh.models.widgets.Slider(
  title='log2 relative transcription rate',
  value=0,
  start=-3,
  end=10,
  step=1)
controls = [rate0_slider, rate1_slider, rate_slider]

prior_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))
obs_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

def update(attr, old, new):
  args = [x.value for x in controls]
  p = np.random.beta(np.exp(np.log(10) * args[0]), np.exp(np.log(10) * args[1]), size=2000)
  g = st.gaussian_kde(p)
  grid = np.linspace(0, 1, 200)
  prior_data.data = bokeh.models.ColumnDataSource.from_df(
    pd.DataFrame({'x': grid, 'y': g(grid)}))

  x = np.random.poisson(np.exp(np.log(2) * args[2]) * p)
  f = st.gaussian_kde(x)
  grid = np.linspace(0, x.max(), 200)
  obs_data.data = bokeh.models.ColumnDataSource.from_df(
    pd.DataFrame({'x': grid, 'y': f(grid)}))

for c in controls:
  c.on_change('value', update)

update(None, None, None)

prior_density = bokeh.plotting.figure(width=300, height=300, tools=[])
prior_density.line(source=prior_data, x='x', y='y', color='black')
prior_density.xaxis.axis_label = 'Prior active proportion'
prior_density.yaxis.axis_label = 'Density'

density = bokeh.plotting.figure(width=300, height=300, tools=[])
density.line(source=obs_data, x='x', y='y', color='black')
density.xaxis.axis_label = 'mRNA count'
density.yaxis.axis_label = 'Density'

layout = bokeh.layouts.layout([[
  bokeh.layouts.widgetbox(*controls, width=300),
  bokeh.layouts.gridplot([[prior_density,  density]])
]])

doc = bokeh.io.curdoc()
doc.title = 'Poisson-Beta simulation'
doc.add_root(layout)
