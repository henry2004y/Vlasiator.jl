# ---
# title: 2D contour plot with streamlines
# id: demo_2d_contour_streamline
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot 2D colored contour with streamlines
# ---

This demos shows how to reading 2D simulation data, zoom-in to a region of interest, and add streamlines on top of colored mesh.
```julia
using Vlasiator, Plots

file = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

boxcoords = Float64[0, 20, -15, 15]

meta = load(file)

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

streamplot(meta, nameV, comp="xy", color="w", density=2.0)
```