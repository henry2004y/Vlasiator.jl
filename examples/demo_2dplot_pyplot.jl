# Sample postprocessing script for customized figure.
#
# Reading 2D simulation data, zooming-in to a region of interest,
# adding streamlines and contour at specific levels on top of colored mesh.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

file = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

meta = load(file)

pcolormesh(meta, nameρ)
streamplot(meta, nameV, comp="xy", color="w", density=2.0)

f, ax = plt.gcf(), plt.gca()
cbar = ax.collections[end].colorbar
boxcoords = [0, 20, -15, 15]
ax.set_xlim([boxcoords[1],boxcoords[2]])
ax.set_ylim([boxcoords[3],boxcoords[4]])

# Contour line at a specific level
pArgs = Vlasiator.set_args(meta, nameρ, RE, Linear)
x, y, data = Vlasiator.plot_prep2d(meta, nameρ, pArgs, "", RE)
CS = plt.contour(x, y, data, levels = [1e7],
                 colors=("k",),linestyles=("-",),linewidths=(0.5,))