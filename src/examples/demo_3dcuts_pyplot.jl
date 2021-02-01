# Sample postprocessing script for customized figure.
#
# Reading 3D simulation data and plotting density in 3 cuts in the same figure.
#
# Hongyang Zhou, hyzhou@umich.edu 02/01/2021

using Vlasiator, PyPlot, PyCall
axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
AxesGrid = axes_grid1.AxesGrid

include("../plot/pyplot.jl")

filename = "bulk1.0001000.vlsv"
nameρ = "proton/vg_rho"

meta = read_meta(filename)

# normal cuts in the x,y,z directions
fig = plt.figure(figsize=(12, 4))
grid = AxesGrid(fig, 111,
                nrows_ncols=(1, 3),
                axes_pad=0.7,
                share_all=false,
                cbar_mode="single",
                cbar_location="right",
                cbar_pad=0.1,
                label_mode="all"
                )

c1 = plot_colormap3dslice(meta, nameρ, grid[1]; normal="x", addcolorbar=false)
c2 = plot_colormap3dslice(meta, nameρ, grid[2]; normal="y", addcolorbar=false)
c3 = plot_colormap3dslice(meta, nameρ, grid[3]; normal="z", addcolorbar=false)

cb = fig.colorbar(c3, cax=grid.cbar_axes[1])
datainfo = read_variable_info(meta, nameρ)

cb_title_str = datainfo.variableLaTeX
cb_title_str *= ",["*datainfo.unitLaTeX*"]"
cb_title = cb.ax.set_title(cb_title_str, fontsize=14, fontweight="bold")

plt.savefig("test.png",bbox_inches="tight")