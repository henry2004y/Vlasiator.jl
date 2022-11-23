# Sample script for mesh plotting.
#
# Represent cell centers as markers.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

const file = "test/data/bulk.amr.vlsv"
meta = load(file)
# 3D mesh
fig = plt.figure()
ax = fig.add_subplot(projection="3d")

plotmesh(meta, marker="+")

# 2D mesh
fig = plt.figure()

pcolormesh(meta, "proton/vg_rho"; axisunit=SI)
plotmesh(meta, projection="y"; color="w")