from netCDF4 import Dataset
import numpy as np

try:
    ncfile.close()
except:
    pass

ncfile = Dataset("test_4.nc", mode="w",
                 format="NETCDF4_CLASSIC")
print(ncfile)

m_dim = 6
n_dim = 4
m = ncfile.createDimension("m", m_dim)
n = ncfile.createDimension("n", n_dim)

for dim in ncfile.dimensions.items():
    print(dim)

ncfile.title='Test for diagonal neighbour'
print(ncfile.title)

mask = ncfile.createVariable('mask', np.int32, ('n', 'm'))

mask[:, :] = np.ones((n_dim, m_dim), dtype=np.int32)