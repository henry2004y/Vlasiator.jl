# ---
# title: 2D contour plots
# id: demo_2d_contour
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to generate 2D colored contours
# ---

This script creates contours from a series of input files with a fixed color range.
Combined with `ffmpeg`, we can easily make animations from data.
```julia
using Vlasiator, Glob, VlasiatorPyPlot, Printf

files = glob("bulk*.vlsv")
nfile = length(files)
# Set output directory
const outdir = "out/"

meta = load(files[1])
const varname = "proton/vg_rho"

fig, ax = plt.subplots()

const comp = :z # vector component for plotting (if applicable)
const axisunit = EARTH
const colorscale = Log
const addcolorbar = true
# Choose colormap
const cmap = matplotlib.cm.turbo
# Set data plotting range
const vmin = 7.0e4
const vmax = 2.5e6

pArgs = Vlasiator.set_args(meta, varname, axisunit; normal=:none)

x1, x2 = Vlasiator.get_axis(pArgs)

norm, ticks = set_colorbar(colorscale, vmin, vmax)

fakedata = zeros(Float32, length(x2), length(x1))

c = ax.pcolormesh(x1, x2, fakedata; norm, cmap)

VlasiatorPyPlot.set_plot(c, ax, pArgs, ticks, addcolorbar)

for (i, file) in enumerate(files)
   @info "$i out of $nfile"
   local meta = load(file)

   var = meta[varname]
   t = readparameter(meta, "time")

   data = Vlasiator.prep2d(meta, varname, comp)'
   c.set_array(data)

   str_title = @sprintf "t= %4.1fs" t
   ax.set_title(str_title)

   savefig(outdir*lpad(i, 4, '0')*".png", bbox_inches="tight", dpi=100)
end

close()
```