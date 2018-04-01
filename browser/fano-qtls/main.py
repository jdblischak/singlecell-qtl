"""Single cell variance QTL browser

"""
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

# This needs to be global to be visible to callbacks
gene = None

def update_gene(attr, old, new):
  selected = gene_data.selected['1d']['indices']
  if not selected:
    return
  with sqlite3.connect(db) as conn:
    global gene
    gene = next(conn.execute('select gene from fano_qtls where fano_qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    print('Selected {}'.format(gene))
    params = pd.read_sql(
      sql="""select fano_qtl_geno.ind, fano_qtl_geno.value as genotype, log_mu,
      log_phi, logodds, mean, var, var / mean as fano from
      fano_qtl_geno, params where fano_qtl_geno.gene == ? and
      fano_qtl_geno.gene == params.gene and fano_qtl_geno.ind == params.ind;""",
      params=(gene,),
      con=conn)
    params['cv'] = np.sqrt(params['var']) / params['mean']
    ind_data.data = bokeh.models.ColumnDataSource.from_df(params)
    umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['left', 'right', 'count']))
    dist_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame(columns=['x', 'y']))

def update_umi(attr, old, new):
  selected = ind_data.selected['1d']['indices']
  with sqlite3.connect(db) as conn:
    if selected:
      ind = ind_data.data['ind'][selected[0]]
      print("Selected {}, {}".format(ind, gene))
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
      sql="""select fano_qtls.gene, gene_info.name, id, p_beta as p, beta from fano_qtls,
      gene_info where fano_qtls.gene == gene_info.gene and fdr_pass order by p_beta;""",
      con=conn))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(
  columns=['ind', 'genotype', 'log_mu', 'log_phi', 'logodds', 'mean', 'var', 'cv', 'fano']))
ind_data.on_change('selected', update_umi)

gene_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['gene', 'id', 'p', 'beta']))
gene_data.on_change('selected', update_gene)

umi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['left', 'right', 'count']))
dist_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

# These need to be module scope because bokeh.server looks there
qtls = bokeh.models.widgets.DataTable(
    source=gene_data,
    columns=[bokeh.models.widgets.TableColumn(field=x, title=x) for x in ['gene', 'name', 'id', 'p', 'beta']],
    width=1200,
    height=300)

hover = bokeh.models.HoverTool(tooltips=[('Individual', '@ind')])

sc_mu_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_mu_by_geno.scatter(source=ind_data, x='genotype', y='log_mu', color='black', size=8)
sc_mu_by_geno.xaxis.axis_label = 'Dosage'
sc_mu_by_geno.yaxis.axis_label = 'log(μ)'

sc_phi_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_phi_by_geno.scatter(source=ind_data, x='genotype', y='log_phi', color='black', size=8)
sc_phi_by_geno.xaxis.axis_label = 'Dosage'
sc_phi_by_geno.yaxis.axis_label = 'log(φ)'

sc_logodds_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_logodds_by_geno.scatter(source=ind_data, x='genotype', y='logodds', color='black', size=8)
sc_logodds_by_geno.xaxis.axis_label = 'Dosage'
sc_logodds_by_geno.yaxis.axis_label = 'logit(π)'

sc_mean_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_mean_by_geno.scatter(source=ind_data, x='genotype', y='mean', color='black', size=8)
sc_mean_by_geno.xaxis.axis_label = 'Dosage'
sc_mean_by_geno.yaxis.axis_label = 'ZI-corrected mean'

sc_var_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_var_by_geno.scatter(source=ind_data, x='genotype', y='var', color='black', size=8)
sc_var_by_geno.xaxis.axis_label = 'Dosage'
sc_var_by_geno.yaxis.axis_label = 'ZI-corrected variance'

sc_cv_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_cv_by_geno.scatter(source=ind_data, x='genotype', y='cv', color='black', size=8)
sc_cv_by_geno.xaxis.axis_label = 'Dosage'
sc_cv_by_geno.yaxis.axis_label = 'ZI-corrected CV'

sc_fano_by_geno = bokeh.plotting.figure(width=300, height=300, tools=['tap', hover])
sc_fano_by_geno.scatter(source=ind_data, x='genotype', y='fano', color='black', size=8)
sc_fano_by_geno.xaxis.axis_label = 'Dosage'
sc_fano_by_geno.yaxis.axis_label = 'ZI-corrected Fano factor'

umi = bokeh.plotting.figure(width=300, height=300, tools=[])
umi.quad(source=umi_data, bottom=0, top='count', left='left', right='right', color='black')
umi.line(source=dist_data, x='x', y='y', color='red', line_width=2)
umi.xaxis.axis_label = 'Observed UMI'
umi.yaxis.axis_label = 'Number of cells'

layout = bokeh.layouts.layout(
  [[qtls],
   [sc_mean_by_geno, sc_var_by_geno, sc_cv_by_geno, sc_fano_by_geno],
   [sc_mu_by_geno, sc_phi_by_geno, sc_logodds_by_geno, umi]],
  sizing_mode='fixed')

doc = bokeh.io.curdoc()
doc.title = 'Fano QTL browser'
doc.add_root(layout)

init()
