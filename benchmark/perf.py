import timeit

p1 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt

file = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
"""

p2 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt

file = 'bulk1.0001000.vlsv'
var = 'proton/vg_rho'
f = pt.vlsvfile.VlsvReader(file)
"""

s1 = """
f = pt.vlsvfile.VlsvReader(file)
"""

s2 = """
cellID = f.read_variable('CellID')
cell_sorted = np.argsort(cellID)
rho = f.read_variable(var)[cell_sorted]
"""

print(timeit.timeit(stmt=s1, setup=p1, number=10) / 10 * 1e3) # [ms]
print(timeit.timeit(stmt=s2, setup=p2, number=10) / 10 * 1e3) # [ms]