import sys
import h5py

data = list()

for l in sys.stdin.readlines():
    data.append(float(l[:-1]))

with h5py.File('test_py.h5', 'w') as f:
    f.create_dataset('distances', data=data)


