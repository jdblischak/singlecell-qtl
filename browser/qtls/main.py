"""Single cell QTL browser

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
    gene = next(conn.execute('select gene from qtls where qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    params = pd.read_sql(
      sql="""select mean_qtl_geno.ind, mean_qtl_geno.value as genotype, log_mu, log_phi,
      logodds, log_mean, log_mean_se, log_phi_se, log_mean + log_mean_se as log_mean_upper, log_mean - log_mean_se as
      log_mean_lower, log_phi + log_phi_se as log_phi_upper, log_phi - log_phi_se
      as log_phi_lower from mean_qtl_geno, params where mean_qtl_geno.gene == ?
      and mean_qtl_geno.gene == params.gene and mean_qtl_geno.ind ==
      params.ind;""",
      params=(gene,),
      con=conn)

    # Compute the regression line
    fit_params = pd.read_sql(sql="""select * from qtls where gene == ?;""", con=conn, params=(gene,))
    x = np.linspace(0, 2, 50)
    y_mean = params['log_mean'].mean() + (x - params['genotype'].mean()) * fit_params.loc[0, 'beta_mean']
    y_disp = params['log_phi'].mean() + (x - params['genotype'].mean()) * fit_params.loc[0, 'beta_disp']
    fit_data.data = bokeh.models.ColumnDataSource.from_df(
      pd.DataFrame({'x': x, 'y_mean': y_mean, 'y_disp': y_disp,
                    'beta_mean': fit_params.loc[0, 'beta_mean'],
                    'pm_mean': fit_params.loc[0, 'pm_mean'],
                    'lfsr_mean': fit_params.loc[0, 'lfsr_mean'],
                    'beta_disp': fit_params.loc[0, 'beta_disp'],
                    'pm_disp': fit_params.loc[0, 'pm_disp'],
                    'lfsr_disp': fit_params.loc[0, 'lfsr_disp'],
      }))

    # Jitter the points
    np.random.seed(0)
    params['genotype'] += np.random.normal(scale=0.01, size=params.shape[0])
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
      sql="""select * from qtls order by lfsr_mean;""",
      con=conn))

# Name this so we can deal with the columns cleanly below
dummy_gene_data = pd.DataFrame(
    columns=['gene', 'name', 'id', 'beta_mean', 'pm_mean', 'lfsr_mean', 'beta_disp', 'pm_disp', 'lfsr_disp',
             'log_mean_resid_var', 'log_mean_error_var', 'log_phi_resid_var', 'log_phi_error_var'])
gene_data = bokeh.models.ColumnDataSource(dummy_gene_data)
gene_data.on_change('selected', update_gene)

fit_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y_mean', 'y_disp']))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(
  columns=['ind', 'genotype', 'log_mu', 'log_phi', 'logodds', 'log_mean', 'log_mean_se', 'log_phi_se', 'log_mean_upper', 'log_mean_lower', 'log_phi_upper', 'log_phi_lower']))
ind_data.on_change('selected', update_umi)

umi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['left', 'right', 'count']))
dist_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

# These need to be module scope because bokeh.server looks there
qtls = bokeh.models.widgets.DataTable(
    source=gene_data,
    columns=[bokeh.models.widgets.TableColumn(field=x, title=x) for x in dummy_gene_data.columns],
    width=1200,
    height=300)

ind_tooltips = [('Individual', '@ind'), ('ln mean SE', '@log_mean_se'), ('ln(φ) SE', '@log_phi_se')]

def generate_fit_tooltips(pheno):
  return [(x, '@{}_{}'.format(x, pheno)) for x in ('beta', 'pm', 'lfsr')]

def generate_figure(pheno, fit_pheno, ylabel):
  res = bokeh.plotting.figure(width=300, height=300, tools=['pan', 'wheel_zoom', 'reset', 'tap'])
  r1 = res.scatter(source=ind_data, x='genotype', y='log_{}'.format(pheno), color='black', size=6)
  r2 = res.segment(source=ind_data, x0='genotype', y0='log_{}_lower'.format(pheno), x1='genotype', y1='log_{}_upper'.format(pheno), color='black', line_width=2)
  res.add_tools(bokeh.models.HoverTool(renderers=[r1, r2], tooltips=ind_tooltips))
  r3 = res.line(source=fit_data, x='x', y='y_{}'.format(fit_pheno), color='red', line_width=2)
  res.add_tools(bokeh.models.HoverTool(renderers=[r3], tooltips=generate_fit_tooltips(fit_pheno)))
  res.xaxis.axis_label = 'Dosage'
  res.yaxis.axis_label = ylabel
  return res
  
sc_log_mean_by_geno = generate_figure(pheno='mean', fit_pheno='mean', ylabel='Deconvolved ln mean')
sc_phi_by_geno = generate_figure(pheno='phi', fit_pheno='disp', ylabel='ln(φ)')

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
