# Sample postprocessing script for customized figure.
#
# This script creates contours from a series of input files with a fixed color range.
# Combined with ffmpeg, we can easily make animations from data.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob, PyPlot, Printf

files = glob("bulk*.vlsv")
nfile = length(files)
# Set output directory
outdir = "out/"

meta = load(files[1])
vardict = Dict("v"=>"proton/vg_v", "rho"=>"proton/vg_rho", "b"=>"vg_b_vol", "b2"=>"fg_b")
varname = "rho"

fig, ax = plt.subplots()

comp = :z # vector component for plotting (if applicable)
axisunit = EARTH
colorscale = Log
addcolorbar = true
# Choose colormap
cmap = matplotlib.cm.turbo
# Set data plotting range
vmin = 7.0e4
vmax = 2.5e6

pArgs = Vlasiator.set_args(meta, vardict[varname], axisunit; normal=:none)

x1, x2 = Vlasiator.get_axis(pArgs)

norm, ticks = Vlasiator.set_colorbar(colorscale, vmin, vmax)

fakedata = zeros(Float32, length(x2), length(x1))

c = ax.pcolormesh(x1, x2, fakedata; norm, cmap)

Vlasiator.set_plot(c, ax, pArgs, ticks, addcolorbar)

for (i, file) in enumerate(files)
   @info "$i out of $nfile"
   local meta = load(file)

   var = meta[vardict[varname]]
   t = readparameter(meta, "time")

   data = Vlasiator.prep2d(meta, vardict[varname], comp)'
   c.set_array(data)

   str_title = @sprintf "t= %4.1fs" t
   ax.set_title(str_title)

   savefig(outdir*lpad(i, 4, '0')*".png", bbox_inches="tight")
end

close()