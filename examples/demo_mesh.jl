# Sample script for mesh plotting.
#
# Represent cell centers as markers.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

filename = "test/data/bulk.amr.vlsv"
meta = load(filename)

fig = plt.figure()

ax = fig.add_subplot(projection="3d")

plotmesh(meta, marker="+")