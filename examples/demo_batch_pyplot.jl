# Sample postprocessing script for customized figure.
#
# Plotting of the same style from a series of input files. Combined with ffmpeg,
# we can easily make animations from data.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob, PyPlot, Printf
using Vlasiator: set_args, plot_prep2d, set_colorbar, set_plot

files = glob("out/bulk*.vlsv")
nfile = length(files)

meta = load(files[1])
vardict = Dict("v"=>"proton/vg_v", "rho"=>"proton/vg_rho", "b"=>"vg_b_vol", "b2"=>"fg_b")
varname = "rho"

fig, ax = plt.subplots()

op = :z
axisunit = RE
colorscale = Log
addcolorbar = true
cmap = matplotlib.cm.turbo
vmin = 7.0e4
vmax = 2.5e6

pArgs = set_args(meta, vardict[varname], axisunit, colorscale; normal=:none, vmin, vmax)

cnorm, cticks = set_colorbar(pArgs)

for (i, file) in enumerate(files)
   @info "$i out of $nfile"
   local meta = load(file)

   var = meta[vardict[varname]]
   t = readparameter(meta, "time")

   if i == 1
      x, y, data = plot_prep2d(meta, vardict[varname], pArgs, op)
   
      c = ax.pcolormesh(x, y, data, norm=cnorm, cmap=cmap, shading="nearest")
   
      set_plot(c, ax, pArgs, cticks, addcolorbar)
   else
      x, y, data = plot_prep2d(meta, vardict[varname], pArgs, op) 
   
      c = ax.pcolormesh(x, y, data, norm=cnorm, cmap=cmap, shading="nearest")
   end

   xlabel("X")
   ylabel("Y")
   str_title = @sprintf "t= %4.1fs" t
   title(str_title)

   savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")

   ax.cla()
end

close()