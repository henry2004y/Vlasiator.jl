# Sample postprocessing script for customized figure.
#
# Reading 2D simulation data, zooming-in to a region of interest,
# adding streamlines and contour at specific levels on top of colored mesh.
#
# Hongyang Zhou, hyzhou@umich.edu 01/25/2021

using Vlasiator, PyPlot

include("../plot/pyplot.jl")

filename = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

meta = read_meta(filename)

plot_pcolormesh(meta, nameρ)
streamline(meta, nameV, comp="xy", color="w", density=2.0)

f, ax = plt.gcf(), plt.gca()
cbar = ax.collections[end].colorbar
boxcoords = [0, 20, -15, 15]
ax.set_xlim([boxcoords[1],boxcoords[2]])
ax.set_ylim([boxcoords[3],boxcoords[4]])

# Contour line at a specific level
pArgs = set_args(meta, nameρ, "Re", true)
x, y, data = plot_prep2d(meta, nameρ, pArgs, "", "Re")
CS = plt.contour(x, y, data, levels = [1e7],
                 colors=("k",),linestyles=("-",),linewidths=(0.5,))