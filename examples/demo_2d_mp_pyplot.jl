# Combined 2D plots across multiple frames, multi-processing version.
#
# To run on a single node,
# julia -p $ncores demo_2d_mp_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, Printf, LaTeXStrings

@everywhere function init_figure(cmaps, norms, ticks, pArgs1, extent)
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

   axs[2,1].set_xlabel(pArgs1.strx; fontsize)
   axs[2,2].set_xlabel(pArgs1.strx; fontsize)
   axs[2,3].set_xlabel(pArgs1.strx; fontsize)
   axs[1,1].set_ylabel(pArgs1.stry; fontsize)
   axs[2,1].set_ylabel(pArgs1.stry; fontsize)

   axs[1,1].set_title("Density"; fontsize)
   axs[1,2].set_title("Velocity x"; fontsize)
   axs[1,3].set_title("Velocity z"; fontsize)
   axs[2,1].set_title("Thermal pressure"; fontsize)
   axs[2,2].set_title("Magnetic field z"; fontsize)
   axs[2,3].set_title("Electric field"; fontsize)

   x1, x2 = Vlasiator.get_axis(pArgs1)

   range1 = searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2])
   range2 = searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])

   fakedata = zeros(Float32, length(range2), length(range1))

   c1 = axs[1,1].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[1], cmap=cmaps[1])
   c2 = axs[1,2].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[2], cmap=cmaps[2])
   c3 = axs[1,3].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[3], cmap=cmaps[2])
   c4 = axs[2,1].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[4], cmap=cmaps[1])
   c5 = axs[2,2].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[5], cmap=cmaps[2])
   c6 = axs[2,3].pcolormesh(x1[range1], x2[range2], fakedata; norm=norms[6], cmap=cmaps[1])

   format = matplotlib.ticker.FormatStrFormatter("%.1f")
   cb1 = colorbar(c1; ax=axs[1,1], ticks=ticks[1], format)
   cb1.ax.set_ylabel("[amu/cc]"; fontsize)
   cb1.outline.set_linewidth(1.0)

   cb2 = colorbar(c2; ax=axs[1,2], ticks=ticks[2])
   cb2.ax.set_ylabel("[km/s]"; fontsize)
   cb2.outline.set_linewidth(1.0)

   cb3 = colorbar(c3; ax=axs[1,3], ticks=ticks[3])
   cb3.ax.set_ylabel("[km/s]"; fontsize)
   cb3.outline.set_linewidth(1.0)

   cb4 = colorbar(c4; ax=axs[2,1], ticks=ticks[4], format)
   cb4.ax.set_ylabel("[nPa]"; fontsize)
   cb4.outline.set_linewidth(1.0)

   cb5 = colorbar(c5; ax=axs[2,2], ticks=ticks[5], format, extend="both")
   cb5.ax.set_ylabel("[nT]"; fontsize)
   cb5.outline.set_linewidth(1.0)

   cb6 = colorbar(c6; ax=axs[2,3], ticks=ticks[6], extend="max")
   cb6.ax.set_ylabel("[mV/m]"; fontsize)
   cb6.outline.set_linewidth(1.0)

   cs = (c1, c2, c3, c4, c5, c6)

   return fig, cs, range1, range2
end

@everywhere function update_plot(fig, cs, range1, range2, file)
   isfile("contour/"*file[end-8:end-5]*".png") && return

   println("file = $file")
   meta = load(file)

   data = Vlasiator.prep2d(meta, "proton/vg_rho", :mag)'
   cs[1].set_array(data[range2,range1] ./ 1f6)

   data = Vlasiator.prep2d(meta, "proton/vg_v", :x)'
   cs[2].set_array(data[range2,range1] ./ 1f3)

   data = Vlasiator.prep2d(meta, "proton/vg_v", :z)'
   cs[3].set_array(data[range2,range1] ./ 1f3)

   data = Vlasiator.prep2d(meta, "vg_pressure", :mag)'
   cs[4].set_array(data[range2,range1] .* 1f9)

   data = Vlasiator.prep2d(meta, "vg_b_vol", :z)'
   cs[5].set_array(data[range2,range1] .* 1f9)

   data = Vlasiator.prep2d(meta, "vg_e_vol", :mag)'
   cs[6].set_array(data[range2,range1] .* 1f3)

   str_title = @sprintf "Density pulse run, t= %4.1fs" meta.time
   fig.suptitle(str_title; fontsize="xx-large")

   savefig("contour/"*file[end-8:end-5]*".png", bbox_inches="tight")
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status, cmaps, norms, ticks, pArgs1, extent)
   fig, cs, range1, range2 = init_figure(cmaps, norms, ticks, pArgs1, extent)
   while true
      file = take!(jobs)
      update_plot(fig, cs, range1, range2, file)
      put!(status, true)
   end
   close(fig)
end

############################################################################################
files = glob("bulk*.vlsv", ".")

nfile = length(files)

# Set colormaps for continuous and divergent data
cmaps = (matplotlib.cm.turbo, matplotlib.cm.RdBu_r)
axisunit = EARTH
extent = [0., 20., -20., 20.] # [RE], default full domain: [-Inf32, Inf32, -Inf32, Inf32]

# Plotting range for each variable
ρmin, ρmax   = 0.1, 15.0     # [amu/cc]
vxmin, vxmax = -650.0, 650.0 # [km/s]
vzmin, vzmax = -500.0, 500.0 # [km/s]
pmin, pmax   = 0.0, 3.6      # [nPa]
bmin, bmax   = -60.0, 60.    # [nT]
emin, emax   = 0.0, 20.      # [mV/m]

meta = load(files[1])

pArgs1 = Vlasiator.set_args(meta, "proton/vg_rho", axisunit; normal=:none)

norm1, ticks1 = Vlasiator.set_colorbar(Linear, ρmin, ρmax)
norm2, ticks2 = Vlasiator.set_colorbar(Linear, vxmin, vxmax)
norm3, ticks3 = Vlasiator.set_colorbar(Linear, vzmin, vzmax)
norm4, ticks4 = Vlasiator.set_colorbar(Linear, pmin, pmax)
norm5, ticks5 = Vlasiator.set_colorbar(Linear, bmin, bmax)
norm6, ticks6 = Vlasiator.set_colorbar(Linear, emin, emax)

norms = (norm1, norm2, norm3, norm4, norm5, norm6)
ticks = (ticks1, ticks2, ticks3, ticks4, ticks5, ticks6)

const fontsize = "x-large"

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

@passobj 1 workers() files

println("Total number of files: $nfile")
println("Running with $(nworkers()) workers...")

@async make_jobs(files) # Feed the jobs channel with all files to process.

@sync for p in workers()
   @async remote_do(do_work, p, jobs, status, cmaps, norms, ticks, pArgs1, extent)
end

let n = nfile
   t = @elapsed while n > 0 # wait for all jobs to complete
      take!(status)
      n -= 1
   end
   println("Finished in $(round(t, digits=2))s.")
end