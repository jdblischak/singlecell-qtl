"""Single cell eQTL browser

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

db = '/scratch/midway2/aksarkar/singlecell/browser.db'

# This needs to be global to be visible to callbacks
gene = None

def update_gene(attr, old, new):
  selected = gene_data.selected['1d']['indices']
  if not selected:
    return
  with sqlite3.connect(db) as conn:
    global gene
    gene = next(conn.execute('select gene from qtls where qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    print('Selected {}'.format(gene))
    ind_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      sql="""select genotype.ind, genotype.value as genotype, params.zinb2_log_mean as
      mean, params.zinb2_log_disp as disp, params.zinb2_logodds as logodds,
      bulk.value as bulk from genotype, bulk, params where genotype.gene == ?
      and genotype.gene == bulk.gene and bulk.gene == params.gene and
      genotype.ind == bulk.ind and bulk.ind == params.ind;""",
      params=(gene,),
      con=conn))

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
      keep = umi['value'] < 19
      edges = np.arange(20)
      counts, _ = np.histogram(umi['value'].values, bins=edges)
      umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame({'left': edges[:-1], 'right': edges[1:], 'count': counts}))

      params = pd.read_sql('select case when nb_nll < zinb_nll then nb_log_mean else zinb2_log_mean end as log_mean, case when nb_nll < zinb_nll then nb_log_disp else zinb2_log_disp end as log_disp, case when nb_nll < zinb_nll then null else zinb2_logodds end as logodds from params where gene == ? and ind == ?', con=conn, params=(gene, ind))
      n = np.exp(params['log_disp'])
      p = 1 / (1 + np.outer(umi['size'], np.exp(params['log_mean'] - params['log_disp'])))
      assert (n > 0).all(), 'n must be non-negative'
      assert (p >= 0).all(), 'p must be non-negative'
      assert (p <= 1).all(), 'p must be <= 1'
      G = st.nbinom(n=n.values.ravel(), p=p.ravel()).pmf
      grid = np.arange(19)
      pmf = np.array([G(x).mean() for x in grid])
      if params['logodds'] is not None:
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
      sql="""select gene_info.gene as gene, gene_info.name, qtls.id, 
      qtls.p_bulk, qtls.beta_bulk, qtls.p_sc, qtls.beta_sc
      from gene_info, qtls 
      where gene_info.gene == qtls.gene
      order by p_bulk;""",
      con=conn))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['ind', 'genotype', 'mean', 'disp', 'logodds', 'bulk']))
ind_data.on_change('selected', update_umi)

gene_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['gene', 'name', 'id', 'p_bulk', 'beta_bulk', 'p_sc', 'beta_sc']))
gene_data.on_change('selected', update_gene)

umi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['left', 'right', 'count']))
dist_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

# These need to be module scope because bokeh.server looks there
qtls = bokeh.models.widgets.DataTable(
    source=gene_data,
    columns=[bokeh.models.widgets.TableColumn(field=x, title=x) for x in ['name', 'id', 'p_bulk', 'beta_bulk', 'p_sc', 'beta_sc']],
    width=1200,
    height=400)

sc_mean_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
sc_mean_by_geno.scatter(source=ind_data, x='genotype', y='mean', color='black', size=8)
sc_mean_by_geno.xaxis.axis_label = 'Centered genotype'
sc_mean_by_geno.yaxis.axis_label = 'Estimated single cell log mean expression'

bulk_mean_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
bulk_mean_by_geno.scatter(source=ind_data, x='genotype', y='bulk', color='black', size=8)
bulk_mean_by_geno.xaxis.axis_label = 'Centered genotype'
bulk_mean_by_geno.yaxis.axis_label = 'Bulk log CPM'

umi = bokeh.plotting.figure(width=400, height=400, tools=[])
umi.quad(source=umi_data, bottom=0, top='count', left='left', right='right', color='black')
umi.line(source=dist_data, x='x', y='y', color='red', line_width=2)
umi.xaxis.axis_label = 'Observed UMI'
umi.yaxis.axis_label = 'Number of cells'

sc_disp_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
sc_disp_by_geno.scatter(source=ind_data, x='genotype', y='disp', color='black', size=8)
sc_disp_by_geno.xaxis.axis_label = 'Centered genotype'
sc_disp_by_geno.yaxis.axis_label = 'Estimated log dispersion'

sc_logodds_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
sc_logodds_by_geno.scatter(source=ind_data, x='genotype', y='logodds', color='black', size=8)
sc_logodds_by_geno.xaxis.axis_label = 'Centered genotype'
sc_logodds_by_geno.yaxis.axis_label = 'Estimated dropout log odds'

layout = bokeh.layouts.layout([[qtls], [bulk_mean_by_geno, sc_mean_by_geno, umi], [sc_disp_by_geno, sc_logodds_by_geno]], sizing_mode='fixed')

doc = bokeh.io.curdoc()
doc.title = 'scQTL browser'
doc.add_root(layout)

init()
