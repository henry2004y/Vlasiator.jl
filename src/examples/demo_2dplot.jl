# Sample postprocessing script for customized figure.
#
# Reading 2D simulation data, zooming-in to a region of interest, and
# adding streamlines on top of colored mesh.
#
# Hongyang Zhou, hyzhou@umich.edu 1/25/2021

using Vlasiator, PyPlot

filename = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

meta = read_meta(filename)

plot_pcolormesh(meta, nameρ)
streamplot(meta, nameV, comp="xy", color="w", density=2.0)

f, ax = plt.gcf(), plt.gca()
boxcoords = [0, 20, -15, 15]
ax.set_xlim([boxcoords[1],boxcoords[2]])
ax.set_ylim([boxcoords[3],boxcoords[4]])
