"""Single cell variance QTL browser

"""
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import numpy as np
import pandas as pd
import scipy.stats as st
import scipy.special as sp
import sqlite3

db = '/project2/mstephens/aksarkar/projects/singlecell-qtl/browser/browser.db'

# This needs to be global to be visible to callbacks
gene = None

def update_gene(attr, old, new):
  selected = gene_data.selected['1d']['indices']
  if not selected:
    return
  with sqlite3.connect(db) as conn:
    global gene
    gene = next(conn.execute('select gene from mean_qtls where mean_qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    params = pd.read_sql(
      sql="""select mean_qtl_geno.ind, mean_qtl_geno.value as genotype, log_mu, log_phi,
      logodds, log_mean, log_mean_se, log_phi_se, log_mean + log_mean_se as log_mean_upper, log_mean - log_mean_se as
      log_mean_lower, log_phi + log_phi_se as log_phi_upper, log_phi - log_phi_se
      as log_phi_lower from mean_qtl_geno, params where mean_qtl_geno.gene == ?
      and mean_qtl_geno.gene == params.gene and mean_qtl_geno.ind ==
      params.ind;""",
      params=(gene,),
      con=conn)
    # Jitter the points
    np.random.seed(0)
    params['genotype'] += np.random.normal(scale=0.01, size=params.shape[0])
    # Scale ln to log2
    for k in params:
      if k.startswith('log'):
        params[k] /= np.log(2)
    ind_data.data = bokeh.models.ColumnDataSource.from_df(params)
    umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['left', 'right', 'count']))
    dist_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['x', 'y']))

def update_umi(attr, old, new):
  selected = ind_data.selected['1d']['indices']
  with sqlite3.connect(db) as conn:
    if selected:
      ind = ind_data.data['ind'][selected[0]]
      umi = pd.read_sql(
        """select umi.value, annotation.size from annotation, umi 
        where umi.gene == ? and annotation.chip_id == ? and 
        umi.sample == annotation.sample""",
        con=conn,
        params=(gene, ind,))
      grid = np.arange(max(umi['value'].values))
      counts, _ = np.histogram(umi['value'].values, bins=grid)
      umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame({'left': grid[:-1], 'right': grid[1:], 'count': counts}))

      params = pd.read_sql('select log_mu, log_phi, logodds from params where gene == ? and ind == ?', con=conn, params=(gene, ind))
      n = np.exp(-params['log_phi'])
      p = 1 / (1 + np.outer(umi['size'], np.exp(params['log_mu'] + params['log_phi'])))
      assert (n > 0).all(), 'n must be non-negative'
      assert (p >= 0).all(), 'p must be non-negative'
      assert (p <= 1).all(), 'p must be <= 1'
      G = st.nbinom(n=n.values.ravel(), p=p.ravel()).pmf
      pmf = np.array([G(x).mean() for x in grid])
      if np.isfinite(params.iloc[0]['logodds']):
        pmf *= sp.expit(-params['logodds']).values
        pmf[0] += sp.expit(params['logodds']).values
      exp_count = umi.shape[0] * pmf
      dist_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame({'x': .5 + grid, 'y': exp_count}))
    else:
      umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['left', 'right', 'count']))
      dist_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['x', 'y']))

def init():
  with sqlite3.connect(db) as conn:
    gene_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      sql="""select gene, name, id, p_beta as p, beta, log_mean_resid_var, log_mean_error_var, log_phi_resid_var, log_phi_error_var from mean_qtls order by p_beta;""",
      con=conn))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(
  columns=['ind', 'genotype', 'log_mu', 'log_phi', 'logodds', 'log_mean', 'log_mean_se', 'log_phi_se', 'log_mean_upper', 'log_mean_lower', 'log_phi_upper', 'log_phi_lower']))
ind_data.on_change('selected', update_umi)

gene_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['gene', 'id', 'p', 'beta', 'log_mean_resid_var', 'log_mean_error_var', 'log_phi_resid_var', 'log_phi_error_var']))
gene_data.on_change('selected', update_gene)

umi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['left', 'right', 'count']))
dist_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

# These need to be module scope because bokeh.server looks there
qtls = bokeh.models.widgets.DataTable(
    source=gene_data,
    columns=[bokeh.models.widgets.TableColumn(field=x, title=x) for x in ['gene', 'name', 'id', 'p', 'beta', 'log_mean_resid_var', 'log_mean_error_var', 'log_phi_resid_var', 'log_phi_error_var']],
    width=1200,
    height=300)

hover = bokeh.models.HoverTool(tooltips=[('Individual', '@ind'), ('log2 mean SE', '@log_mean_se'), ('log2(φ) SE', '@log_phi_se')])

sc_phi_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['pan', 'wheel_zoom', 'reset', 'tap', hover])
sc_phi_by_geno.scatter(source=ind_data, x='genotype', y='log_phi', color='black', size=6)
sc_phi_by_geno.segment(source=ind_data, x0='genotype', y0='log_phi_lower', x1='genotype', y1='log_phi_upper', color='black', line_width=2)
sc_phi_by_geno.xaxis.axis_label = 'Dosage'
sc_phi_by_geno.yaxis.axis_label = 'log2(φ)'

sc_log_mean_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['pan', 'wheel_zoom', 'reset', 'tap', hover])
sc_log_mean_by_geno.scatter(source=ind_data, x='genotype', y='log_mean', color='black', size=6)
sc_log_mean_by_geno.segment(source=ind_data, x0='genotype', y0='log_mean_lower', x1='genotype', y1='log_mean_upper', color='black', line_width=2)
sc_log_mean_by_geno.xaxis.axis_label = 'Dosage'
sc_log_mean_by_geno.yaxis.axis_label = 'Deconvolved log2 mean'

umi = bokeh.plotting.figure(width=300, height=300, tools=[])
umi.quad(source=umi_data, bottom=0, top='count', left='left', right='right', color='black')
umi.line(source=dist_data, x='x', y='y', color='red', line_width=2)
umi.xaxis.axis_label = 'Observed UMI'
umi.yaxis.axis_label = 'Number of cells'

panels = bokeh.layouts.gridplot([
   [sc_log_mean_by_geno, sc_phi_by_geno, umi],
])

layout = bokeh.layouts.layout([[qtls], [panels]], sizing_mode='fixed')

doc = bokeh.io.curdoc()
doc.title = 'scQTL browser'
doc.add_root(layout)

init()
