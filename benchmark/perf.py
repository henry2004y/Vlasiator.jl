import timeit

p0 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
"""

p1 = p0+"""
file = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
"""

p2 = p0+"""
file = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
f = pt.vlsvfile.VlsvReader(file)
"""
# Reading meta data
s1 = """
f = pt.vlsvfile.VlsvReader(file)
"""
# Obtaining scalar variable
s2 = """
cellID = f.read_variable('CellID')
cell_sorted = np.argsort(cellID)
rho = f.read_variable(var)[cell_sorted]
"""
# Plotting contour from uniform mesh
s3 = """
pt.plot.plot_colormap(filename='bulk.0000501.vlsv', var='rho', draw=1)
"""
# Plotting contour slice from AMR mesh
s4 = """
pt.plot.plot_colormap3dslice(filename='bulk1.0001000.vlsv', var='proton/vg_rho', draw=1, normal='y')
"""

p3 = p0+"""
import glob
dir = "/wrk/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"
var = "proton/vg_rho"
Re = 6.371e6 # Earth radius, [m]
loc = [12*Re, 0, 0]
filenames = sorted(glob.glob(dir+"bulk1*.vlsv"))
data_series = np.zeros(len(filenames))
f = pt.vlsvfile.VlsvReader(filenames[0])
id = f.get_cellid(loc)
print(id)
"""
# Extracting variable from static location
s5 = """
for (i,fname) in enumerate(filenames[0:200]):
   f = pt.vlsvfile.VlsvReader(fname)
   data_series[i] = f.read_variable(name=var, cellids=id)
"""

print(timeit.timeit(stmt=s1, setup=p1, number=10) / 10 * 1e3) # [ms]
print(timeit.timeit(stmt=s2, setup=p2, number=10) / 10 * 1e3) # [ms]

print(timeit.timeit(stmt=s3, setup=p0, number=10) / 10) # [s]
print(timeit.timeit(stmt=s4, setup=p0, number=10) / 10) # [s]

print(timeit.timeit(stmt=s5, setup=p3, number=10) / 10) # [s]