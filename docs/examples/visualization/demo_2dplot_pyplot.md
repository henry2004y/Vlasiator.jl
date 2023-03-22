# ---
# title: 2D contour plot with streamlines and levels
# id: demo_2d_contour_streamline_levels
# cover: ../../src/figures/magnetosphere_earth_proton_density_2D.png
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot 2D colored contour with streamlines and levels
# ---

This demo shows how to plot 2D colored contour in a region of interest with streamlines and levels.
```julia
using Vlasiator, PyPlot, VlasiatorPyPlot
using Vlasiator: RE

function main(file::String, varname::String)
   meta = load(file)

   fig, ax = plt.subplots(1,1; figsize=(8,6), constrained_layout=true)

   pcolormesh(meta, varname, ax)
   streamplot(meta, "rho_v", ax; comp="xy", color="w", density=2.0)

   cbar = ax.collections[end].colorbar
   boxcoords = Float64[0, 20, -15, 15]
   ax.set_xlim([boxcoords[1],boxcoords[2]])
   ax.set_ylim([boxcoords[3],boxcoords[4]])

   # Contour line at a specific level
   pArgs = Vlasiator.set_args(meta, varname, EARTH)
   x, y = Vlasiator.get_axis(pArgs)
   data = Vlasiator.prep2d(meta, varname)'
   CS = plt.contour(x, y, data, levels = [1e7],
                    colors=("k",),linestyles=("-",),linewidths=(0.5,))

   # Add a rectangular box region
   boxrange = (250:299, 200:249)

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
   y = LinRange(meta.coordmin[2], meta.coordmax[2], meta.ncells[2])
   # Construct the patch
   rect = matplotlib.patches.Rectangle((x[boxrange[1][1]]/RE, y[boxrange[2][1]]/RE),
      (x[boxrange[1][end]] - x[boxrange[1][1]])/RE,
      (y[boxrange[2][end]] - y[boxrange[2][1]])/RE,
      linewidth=1, edgecolor="r", facecolor="none")

   ax.add_patch(rect)
end


file = "bulk.0000501.vlsv"
nameρ = "rho"

main(file, nameρ)
```