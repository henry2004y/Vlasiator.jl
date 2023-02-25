# ---
# title: Plot mesh
# id: demo_plot_mesh
# cover: ../../src/figures/mesh.png
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot mesh
# ---

This demo shows how to plot mesh. The cell centers are represented by the markers.
```julia
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
```