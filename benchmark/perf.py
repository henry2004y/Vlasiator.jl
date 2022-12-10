# Script for generating benchmark results in
# https://henry2004y.github.io/Vlasiator.jl/dev/benchmark/
#
# Note: this script is not intended for execution as a whole: each benchmark
# needs to be executed independently after importing all the required packages.

import timeit

# Setup: importing packages
p0 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
"""
# Setup: selecting file and variable
p1 = p0+"""
file1 = 'bulk.singleprecision.vlsv' # 80 B
file2 = 'bulk.0000003.vlsv'         # 900 KB
file3 = 'bulk1.0001000.vlsv'        # 32 MB

file = file1
var = 'proton/vg_rho'
"""
# Setup: reading metadata
p2 = p0+p1+"""
f = pt.vlsvfile.VlsvReader(file)
"""
# Setup: virtual satellite extraction
p3 = p0+"""
import glob
dir = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"
var = "proton/vg_rho"
Re = 6.371e6 # Earth radius, [m]
loc = [12*Re, 0, 0]
filenames = sorted(glob.glob(dir+"bulk1*.vlsv"))
data_series = np.zeros(len(filenames))
f = pt.vlsvfile.VlsvReader(filenames[0])
id = f.get_cellid(loc)
print(id)
"""

# Run: reading metadata
s1 = """
f = pt.vlsvfile.VlsvReader(file)
"""
# Run: obtaining scalar variable
s2 = """
cellID = f.read_variable('CellID')
cell_sorted = np.argsort(cellID)
rho = f.read_variable(var)[cell_sorted]
"""
# Run: plotting density contour on a uniform mesh
s3 = """
pt.plot.plot_colormap(filename='bulk.0000501.vlsv', var='rho', draw=1)
"""
# Run: plotting density slice from a 3D AMR mesh
s4 = """
pt.plot.plot_colormap3dslice(filename='bulk1.0001000.vlsv', var='proton/vg_rho', draw=1, normal='y')
"""
# Run: virtual satellite extraction from a static location
s5 = """
for (i,fname) in enumerate(filenames):
   f = pt.vlsvfile.VlsvReader(fname)
   data_series[i] = f.read_variable(name=var, cellids=id)
"""

print(timeit.timeit(stmt=s1, setup=p1, number=10) / 10 * 1e3) # [ms]
print(timeit.timeit(stmt=s2, setup=p2, number=10) / 10 * 1e3) # [ms]

t1 = timeit.timeit(stmt=s1, setup=p1, number=10) / 10
print(f"Finished reading metadata in {t1:0.4f} ms")

t2 = timeit.timeit(stmt=s2, setup=p2, number=10) / 10
print(f"Finished reading scalar DCCRG variable in {t2:0.4f} ms")

t3 = timeit.timeit(stmt=s3, setup=p0, number=3) / 3
print(f"Finished plotting in {t3:0.4f} s")

t4 = timeit.timeit(stmt=s4, setup=p0, number=3) / 3
print(f"Finished plotting in {t4:0.4f} s")

t5 = timeit.timeit(stmt=s5, setup=p3, number=1) / 1
print(f"Finished virtual satellite tracking in {t5:0.4f} s")