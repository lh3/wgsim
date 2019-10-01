import sys
import h5py
import pandas as pd
from scipy.spatial.distance import squareform

distmat = None
names = None
with h5py.File(sys.argv[1], 'r') as f:
    distmat = f['distances'][:]
    names = f['leaf_names'][:].astype('U')

print(distmat)
print(names)
sq = squareform(distmat)

df = pd.DataFrame(data=sq, columns=names, index=names)

print(df)
