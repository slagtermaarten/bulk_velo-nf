import os
import numpy as np
import pandas as pd
import scvelo as scv
import sys

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

fn = sys.argv[1]
adata = scv.read(fn, cache=True)

for t in ['ambiguous', 'matrix', 'spliced', 'unspliced']:
  M = adata.layers[t].toarray().transpose()
  df = pd.DataFrame(M)
  df.index = adata.var.index
  df.columns = adata.obs.index.values
  df.to_csv(os.path.join(os.path.dirname(fn), '{}.csv'.format(t)))
