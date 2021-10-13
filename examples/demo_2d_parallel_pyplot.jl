# Combined 2D plots across multiple frames, multi-process version.
#
# To run on a single node,
# julia -p $ncores demo_2d_parallel_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, Printf, LaTeXStrings
@everywhere using Vlasiator: set_args, plot_prep2d, set_colorbar, set_plot

@everywhere function init_figure()
   fig, axs = plt.subplots(2, 3; num=myid(),
      figsize=(13.5, 9.5), sharex=true, sharey=true, constrained_layout=true)

   for ax in axs
      ax.set_aspect("equal")
      # Set border line widths
      for loc in ("left", "bottom", "right", "top")
         edge = get(ax.spines, loc, nothing)
         edge.set_linewidth(2.0)
      end
      ax.xaxis.set_tick_params(width=2.0, length=3)
      ax.yaxis.set_tick_params(width=2.0, length=3)
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   end

   axs[2,1].set_xlabel(L"X [$R_E$]"; fontsize)
   axs[2,2].set_xlabel(L"X [$R_E$]"; fontsize)
   axs[2,3].set_xlabel(L"X [$R_E$]"; fontsize)
   axs[1,1].set_ylabel(L"Y [$R_E$]"; fontsize)
   axs[2,1].set_ylabel(L"Y [$R_E$]"; fontsize)

   axs[1,1].set_title("Density"; fontsize)
   axs[1,2].set_title("Velocity x"; fontsize)
   axs[1,3].set_title("Velocity y"; fontsize)
   axs[2,1].set_title("Thermal pressure"; fontsize)
   axs[2,2].set_title("Magnetic field"; fontsize)
   axs[2,3].set_title("Electric field"; fontsize)

   plotrange, sizes, axisunit = pArgs1.plotrange, pArgs1.sizes, pArgs1.axisunit
   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   fakedata = zeros(Float32, length(y), length(x))
   c1 = axs[1,1].pcolormesh(x, y, fakedata; norm=cnorm1, cmap=cmap, shading="nearest")
   c2 = axs[1,2].pcolormesh(x, y, fakedata; norm=cnorm2, cmap=cmap, shading="nearest")
   c3 = axs[1,3].pcolormesh(x, y, fakedata;
      norm=cnorm3, cmap=matplotlib.cm.RdBu, shading="nearest")
   c4 = axs[2,1].pcolormesh(x, y, fakedata; norm=cnorm4, cmap=cmap, shading="nearest")
   c5 = axs[2,2].pcolormesh(x, y, fakedata; norm=cnorm5, cmap=cmap, shading="nearest")
   c6 = axs[2,3].pcolormesh(x, y, fakedata; norm=cnorm6, cmap=cmap, shading="nearest")

   format = matplotlib.ticker.FormatStrFormatter("%.1f")
   cb1 = colorbar(c1; ax=axs[1,1], ticks=cticks1, format)
   cb1.ax.set_ylabel("[amu/cc]"; fontsize)
   cb1.outline.set_linewidth(1.0)

   cb2 = colorbar(c2; ax=axs[1,2], ticks=cticks2)
   cb2.ax.set_ylabel("[km/s]"; fontsize)
   cb2.outline.set_linewidth(1.0)

   cb3 = colorbar(c3; ax=axs[1,3], ticks=cticks3)
   cb3.ax.set_ylabel("[km/s]"; fontsize)
   cb3.outline.set_linewidth(1.0)

   cb4 = colorbar(c4; ax=axs[2,1], ticks=cticks4, format)
   cb4.ax.set_ylabel("[nPa]"; fontsize)
   cb4.outline.set_linewidth(1.0)

   cb5 = colorbar(c5; ax=axs[2,2], ticks=cticks5, format)
   cb5.ax.set_ylabel("[nT]"; fontsize)
   cb5.outline.set_linewidth(1.0)

   cb6 = colorbar(c6; ax=axs[2,3], ticks=cticks6)
   cb6.ax.set_ylabel(L"[$\mu V/m$]"; fontsize)
   cb6.outline.set_linewidth(1.0)

   cs = (c1, c2, c3, c4, c5, c6)

   return fig, cs
end

@everywhere function process(fig, cs, file)
   isfile("out/"*file[end-8:end-5]*".png") && return

   println("file = $file")
   meta = load(file)

   _, _, data = plot_prep2d(meta, "proton/vg_rho", pArgs1, :mag)
   cs[1].set_array(data ./ 1e6)

   _, _, data = plot_prep2d(meta, "proton/vg_v", pArgs2, :x)
   cs[2].set_array(data ./ 1e3)

   _, _, data = plot_prep2d(meta, "proton/vg_v", pArgs3, :y)
   cs[3].set_array(data ./ 1e3)

   _, _, data = plot_prep2d(meta, "vg_pressure", pArgs4, :mag)
   cs[4].set_array(data .* 1e9)

   _, _, data = plot_prep2d(meta, "vg_b_vol", pArgs5, :z)
   cs[5].set_array(data .* 1e9)

   _, _, data = plot_prep2d(meta, "vg_e_vol", pArgs6, :mag)
   cs[6].set_array(data .* 1e6)

   str_title = @sprintf "Density pulse run, t= %4.1fs" meta.time
   fig.suptitle(str_title; fontsize="xx-large")

   savefig("out/"*file[end-8:end-5]*".png", bbox_inches="tight")
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status)
   fig, cs = init_figure()
   while true
      file = take!(jobs)
      process(fig, cs, file)
      put!(status, true)
   end
   close(fig)
end

############################################################################################
files = glob("bulk*.vlsv", ".")

nfile = length(files)

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

@passobj 1 workers() files

@broadcast begin # on all workers
   # Set contour plots' axes and colorbars
   const cmap = matplotlib.cm.turbo
   colorscale = Linear
   axisunit = RE

   # Upper/lower limits for each variable
   const ρmin, ρmax   = 0.0, 9.0      # [amu/cc]
   const vxmin, vxmax = -640.0, 100.0 # [km/s]
   const vymin, vymax = -300.0, 300.0 # [km/s]
   const pmin, pmax   = 0.0, 0.7      # [nPa]
   const bmin, bmax   = -5.0, 50.     # [nT]
   const emin, emax   = 5.0, 3e4      # [muV/m]

   meta = load(files[1])

   const pArgs1 = set_args(meta, "proton/vg_rho", axisunit, colorscale;
      normal=:none, vmin=ρmin, vmax=ρmax)
   const cnorm1, cticks1 = set_colorbar(pArgs1)

   const pArgs2 = set_args(meta, "proton/vg_v", axisunit, colorscale;
      normal=:none, vmin=vxmin, vmax=vxmax)
   const cnorm2, cticks2 = set_colorbar(pArgs2)

   const pArgs3 = set_args(meta, "proton/vg_v", axisunit, colorscale;
      normal=:none, vmin=vymin, vmax=vymax)
   const cnorm3, cticks3 = set_colorbar(pArgs3)

   const pArgs4 = set_args(meta, "vg_pressure", axisunit, colorscale;
      normal=:none, vmin=pmin, vmax=pmax)
   const cnorm4, cticks4 = set_colorbar(pArgs4)

   const pArgs5 = set_args(meta, "vg_b_vol", axisunit, Linear;
      normal=:none, vmin=bmin, vmax=bmax)
   const cnorm5, cticks5 = set_colorbar(pArgs5)

   const pArgs6 = set_args(meta, "vg_e_vol", axisunit, Log;
      normal=:none, vmin=emin, vmax=emax)
   const cnorm6, cticks6 = set_colorbar(pArgs6)

   const fontsize = "x-large"
end

println("Total number of files: $nfile")
println("Running with $(nworkers()) workers...")

@async make_jobs(files) # Feed the jobs channel with all files to process.

@sync for p in workers()
   @async remote_do(do_work, p, jobs, status)
end

let n = nfile
   t = @elapsed while n > 0 # wait for all jobs to complete
      take!(status)
      n -= 1
   end
   println("Finished in $(round(t, digits=2))s.")
end