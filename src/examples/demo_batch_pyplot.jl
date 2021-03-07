# Sample postprocessing script for customized figure.
#
# Plotting of the same style from a series of input files. Combined with ffmpeg,
# we can easily make animations from data.
#
# Hongyang Zhou, hyzhou@umich.edu 03/07/2021

using Vlasiator, Glob, PyPlot, Printf

filenames = glob("BCtest/run8_coarse/bulk*.vlsv")
nfile = length(filenames)

meta = read_meta(filenames[1])
varnames = show_variables(meta)
vardict = Dict("v"=>"proton/vg_v", "rho"=>"proton/vg_rho", "b"=>"vg_b_vol",
   "b2"=>"fg_b")
varname = "rho"

fig, ax = plt.subplots()

op= :z
axisunit=""
islinear = false
addcolorbar = true
cmap = matplotlib.cm.turbo
vmin = 7.0e4
vmax = 2.5e6

if islinear
   nticks = 7
   levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
   cnorm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=true)
   cticks = range(vmin, vmax, length=nticks)
else
   cnorm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
   cticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
end

pArgs = set_args(meta, vardict[varname], axisunit, islinear)

for (i, filename) in enumerate(filenames)
   @info "$i out of $nfile"
   local meta = read_meta(filename)

   var = read_variable(meta, vardict[varname])
   t = read_parameter(meta, "time")

   if i == 1
      x, y, data = plot_prep2d(meta, vardict[varname], pArgs, op, axisunit)
   
      c = ax.pcolormesh(x, y, data, norm=cnorm, cmap=cmap, shading="auto")
   
      set_plot(c, ax, pArgs, cticks, addcolorbar)
   else
      x, y, data = plot_prep2d(meta, vardict[varname], pArgs, op, axisunit) 
   
      c = ax.pcolormesh(x, y, data, norm=cnorm, cmap=cmap, shading="auto")
   end

   xlabel("X")
   ylabel("Y")
   str_title = @sprintf "t= %4.1fs" t
   title(str_title)

   savefig("BCtest/"*lpad(i, 4, '0')*".png", bbox_inches="tight")

   ax.cla()
end

close()