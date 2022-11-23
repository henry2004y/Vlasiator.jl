# Sample postprocessing script for customized figure.
#
# Reading 2D simulation data, zooming-in to a region of interest,
# adding streamlines and contour at specific levels on top of colored mesh.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot
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