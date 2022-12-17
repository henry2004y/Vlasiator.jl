# Script for generating benchmark results in
# https://henry2004y.github.io/Vlasiator.jl/dev/benchmark/
#
# If Tex is not installed, LaTeX format output can be disabled via
#    export PTNOLATEX=1

# Set environment variables for Analysator
import os

os.environ["PTNONINTERACTIVE"] = "1"
os.environ["PTBACKEND"] = 'Agg'

import requests, tarfile, timeit, os.path, textwrap
import pytools as pt

files = ['1d_single.vlsv', 'bulk.2d.vlsv', '2d_double.vlsv', '2d_AFC.vlsv', '3d_EGI.vlsv']

for i, file in enumerate(files):
   if os.path.isfile(file):
      print(f"Benchmark file {file} found...")
   elif i == 3:
      url_base = 'https://a3s.fi/swift/v1/AUTH_81f1cd490d494224880ea77e4f98490d/vlasiator-2d-afc/'
      filename = 'production_halfres/bulk.0000000.vlsv'
      r = requests.get(url_base+filename, allow_redirects=True)
      open('2d_AFC.vlsv', 'wb').write(r.content)
   elif i in [0,1,2]:
      r = requests.get('https://raw.githubusercontent.com/henry2004y/vlsv_data/master/testdata.tar.gz', allow_redirects=True)
      open('testdata.tar.gz', 'wb').write(r.content)
      r = requests.get('https://raw.githubusercontent.com/henry2004y/vlsv_data/master/1d_single.vlsv', allow_redirects=True)
      open('1d_single.vlsv', 'wb').write(r.content)
      r = requests.get('https://raw.githubusercontent.com/henry2004y/vlsv_data/master/2d_double.vlsv', allow_redirects=True)
      open('2d_double.vlsv', 'wb').write(r.content)

      file = tarfile.open('testdata.tar.gz')
      file.extractall('./')
      file.close()
   elif i == 4:
      print(f"{file} is not open-access!")

# Number of trails for each benchmark
ntrail = 100

# Setup: importing packages
p0 = """
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
"""

for i, file in enumerate(files):
   if not os.path.isfile(file):
      continue
   # Setup: reading metadata
   p2 = textwrap.dedent("""\
   f = pt.vlsvfile.VlsvReader(file)
   cellID = f.read_variable('CellID')
   cell_sorted = np.argsort(cellID)
   """)

   if i != 3:
      # Setup: selecting file and variable
      p1 = p0 + textwrap.dedent(f"""\
      file = '{file}'
      var = 'proton/vg_rho'
      """)
   else:
      # Setup: selecting file and variable
      p1 = p0 + textwrap.dedent(f"""\
      file = '{file}'
      var = 'proton/rho'
      """)

   # Run: reading metadata
   s1 = p2
   t1 = timeit.timeit(stmt=s1, setup=p1, number=ntrail) / ntrail * 1e3
   print(f"{file}:")
   print(f"Reading metadata in {t1:0.4f} ms")
   # Run: obtaining sorted scalar variable
   s2 = textwrap.dedent("""\
   rho = f.read_variable(var)[cell_sorted]
   """)
   t2 = timeit.timeit(stmt=s2, setup=p0+p1+p2, number=ntrail) / ntrail * 1e3
   print(f"{file}:")
   print(f"Reading scalar DCCRG variable in {t2:0.4f} ms")

# Run: plotting density contour on a uniform mesh
s3 = f"""
pt.plot.plot_colormap(filename='{files[3]}', var='rho', draw=1)
"""
t3 = timeit.timeit(stmt=s3, setup=p0, number=5) / 5
print(f"Uniform 2d plotting in {t3:0.4f} s")

if os.path.isfile(files[4]):
   # Run: plotting density slice from a 3D AMR mesh
   s4 = textwrap.dedent(f"""\
   pt.plot.plot_colormap3dslice(filename='{files[4]}', var='proton/vg_rho', draw=1, normal='y')
   """)

   t4 = timeit.timeit(stmt=s4, setup=p0, number=5) / 5
   print(f"AMR slice plotting in {t4:0.4f} s")