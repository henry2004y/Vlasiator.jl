# Sample script for mesh plotting.
#
# Represent cell centers as markers.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

const file = "test/data/bulk.amr.vlsv"
meta = load(file)

fig = plt.figure()

ax = fig.add_subplot(projection="3d")

plotmesh(meta, marker="+")