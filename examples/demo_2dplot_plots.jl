# Sample postprocessing script for customized figure.
#
# Reading 2D simulation data, zooming-in to a region of interest, and
# adding streamlines on top of colored mesh.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Plots

filename = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

boxcoords = [0, 20, -15, 15]

meta = readmeta(filename)

heatmap(meta, nameρ, 
   xlim=(boxcoords[1], boxcoords[2]),
   ylim=(boxcoords[3], boxcoords[4]),
   aspect_ratio=:equal,
   c=:turbo)

#=
# Attributes can be modified afterwards, but it's slower.
heatmap(meta, nameρ, c=:turbo)
xlims!(boxcoords[1], boxcoords[2])
ylims!(boxcoords[3], boxcoords[4])
=#

#streamplot(meta, nameV, comp="xy", color="w", density=2.0)
