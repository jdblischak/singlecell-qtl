import numpy as np
import pandas as pd
import sqlite3
import os.path

outfile = os.path.join(os.path.dirname(__file__), 'browser.db')

gene_info = (pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-genes.txt.gz')
             .set_index('gene')
             .query('source == "H. sapiens"')
             .query('chr != "hsX"')
             .query('chr != "hsY"')
             .query('chr != "hsMT"'))
with sqlite3.connect(outfile) as conn:
  gene_info.to_sql(name='gene_info', con=conn, if_exists='replace')

def melt_write(df, conn, **kwargs):
  default_kwargs = {'index': False, 'if_exists': 'replace'}
  default_kwargs.update(kwargs)
  (df
   .reset_index()
   .melt(id_vars='gene', var_name='ind')
   .to_sql(con=conn, **default_kwargs))

mean = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean.txt.gz', index_col='gene', sep=' ')
disp = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-dispersion.txt.gz', index_col='gene', sep=' ')
dropout = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-dropout.txt.gz', index_col='gene', sep=' ')

with sqlite3.connect(outfile) as conn:
  melt_write(mean, conn, name='log_mean')
  melt_write(disp, conn, name='log_disp')
  melt_write(dropout, conn, name='logodds')

genotypes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/bulk-qtl-genotypes.txt.gz', index_col='gene', sep=' ')
with sqlite3.connect(outfile) as conn:
  (genotypes
   .reset_index()
   .melt(id_vars='gene', var_name='ind').to_sql(name='genotype', con=conn, index=False, if_exists='replace'))

bulk = (pd.read_table('/project2/gilad/singlecell-qtl/bulk/counts_RNAseq_iPSC.txt', sep=' ', index_col='gene')
        .rename(index=lambda x: x.split('.')[0], columns=lambda x: 'NA{}'.format(x))
        .align(mean, axis=None, join='inner')[0])
bulk = np.log((bulk + 1) / bulk.sum(axis=0))
with sqlite3.connect(outfile) as conn:
  melt_write(bulk, conn, name='bulk')

with sqlite3.connect(outfile) as conn:
  (pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean-qtls.txt.gz', sep=' ', index_col='gene')
   .to_sql('qtls', con=conn, if_exists='replace'))

annotations = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-annotation.txt')
keep_genes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/genes-pass-filter.txt', index_col=0, header=None)
keep_samples = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/quality-single-cells.txt', index_col=0, header=None)
annotations = annotations.loc[keep_samples.values.ravel()]
annotations['sample'] = annotations.apply(lambda x: '.'.join([x['chip_id'], '{:08d}'.format(x['experiment']), x['well']]), axis=1)
annotations = annotations.set_index('sample')
annotations['size'] = np.zeros(annotations.shape[0])
with sqlite3.connect(outfile) as conn:
  for chunk in pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-counts.txt.gz', index_col=0, chunksize=100):
    chunk = (chunk
             .filter(items=keep_genes[keep_genes.values.ravel()].index, axis='index')
             .loc[:,keep_samples.values.ravel()])
    if not chunk.empty:
      annotations['size'] += chunk.sum(axis=0)
      chunk.reset_index().melt(id_vars='gene', var_name='sample').to_sql(name='umi', con=conn, index=False, if_exists='append')
  annotations[['chip_id', 'size']].to_sql(name='annotation', con=conn, if_exists='replace')
