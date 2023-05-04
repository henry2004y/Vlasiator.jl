# ---
# title: Orthogonal slices
# id: demo_3d_cuts
# cover: ../../src/figures/magnetosphere_earth_proton_density_3cuts.png
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot 2D colored contour slices for 3D data
# ---

This demo shows how to plot 2D colored contour slices for 3D data.
```julia
using VlasiatorPyPlot, PyCall

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
ImageGrid = axes_grid1.ImageGrid

file = "bulk1.0001000.vlsv"
nameρ = "proton/vg_rho"
colorscale = Log
addcolorbar = false

meta = load(file)

# normal cuts in the x,y,z directions
fig = plt.figure(figsize=(12, 4))
grid = ImageGrid(fig, 111,
                 nrows_ncols=(1, 3),
                 axes_pad=0.85,
                 cbar_mode="single",
                 cbar_location="right",
                 cbar_pad=0.1,
                 label_mode="all"
                 )

c1 = pcolormesh(meta, nameρ, grid[1]; normal=:x, addcolorbar, colorscale)
c2 = pcolormesh(meta, nameρ, grid[2]; normal=:y, addcolorbar, colorscale)
c3 = pcolormesh(meta, nameρ, grid[3]; normal=:z, addcolorbar, colorscale)

cb = fig.colorbar(c3, cax=grid.cbar_axes[1])
datainfo = readvariablemeta(meta, nameρ)

cb_title_str = datainfo.variableLaTeX
cb_title_str *= ",["*datainfo.unitLaTeX*"]"
cb_title = cb.ax.set_title(cb_title_str, fontsize=14, fontweight="bold")

plt.savefig("test.png", bbox_inches="tight")
```