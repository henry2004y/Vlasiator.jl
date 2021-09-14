import timeit

p1 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt

filename = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
"""

p2 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt

filename = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
f = pt.vlsvfile.VlsvReader(filename)
"""

s1 = """
f = pt.vlsvfile.VlsvReader(filename)
"""

s2 = """
cellID = f.read_variable('CellID')
cell_sorted = np.argsort(cellID)
rho = f.read_variable(var)[cell_sorted]
"""

print(timeit.timeit(stmt=s1, setup=p1, number=10) / 10 * 1e6) # [μs]
print(timeit.timeit(stmt=s2, setup=p2, number=10) / 10 * 1e6) # [μs]