# ---
# title: 2D contour plots with multi-processing
# id: demo_2d_contour_mp
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to generate 2D colored contours with multi-processing
# ---

This demo shows how to generate 2D colored contours with multi-processing.
To run on a single node,
```shell
julia -p $ncores demo_2d_mp_pyplot.jl
```

```julia
using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, VlasiatorPyPlot, Printf, LaTeXStrings

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

   for ax in axs[end,:]
      ax.set_xlabel(pArgs1.strx; fontsize)
   end
   for ax in axs[:,1]
      ax.set_ylabel(pArgs1.stry; fontsize)
   end

   titles = ("Density", "Pth", "Vx", "Bz", "Vz", "E")

   for (ax, title) in zip(axs, titles)
      ax.set_title(title; fontsize)
   end

   x1, x2 = Vlasiator.get_axis(pArgs1)

   range1 = searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2])
   range2 = searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])

   fakedata = zeros(Float32, length(range2), length(range1))

   c1 = axs[1,1].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[1], cmap=cmaps[1])
   c2 = axs[1,2].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[2], cmap=cmaps[2])
   c3 = axs[1,3].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[3], cmap=cmaps[2])
   c4 = axs[2,1].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[4], cmap=cmaps[1])
   c5 = axs[2,2].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[5], cmap=cmaps[2])
   c6 = axs[2,3].pcolormesh(x1[range1], x2[range2], fakedata;
      norm=norms[6], cmap=cmaps[1])

   cs = (c1, c2, c3, c4, c5, c6)

   format = matplotlib.ticker.FormatStrFormatter("%.1f")

   cb1 = colorbar(c1; ax=axs[1,1], ticks=ticks[1], format)
   cb2 = colorbar(c2; ax=axs[1,2], ticks=ticks[2])
   cb3 = colorbar(c3; ax=axs[1,3], ticks=ticks[3])
   cb4 = colorbar(c4; ax=axs[2,1], ticks=ticks[4], format)
   cb5 = colorbar(c5; ax=axs[2,2], ticks=ticks[5], format, extend="both")
   cb6 = colorbar(c6; ax=axs[2,3], ticks=ticks[6], extend="max")

   ylabels = ("[amu/cc]", "[km/s]", "[km/s]", "[nPa]", "[nT]", "[mV/m]")

   for (i, cb) in enumerate((cb1, cb2, cb3, cb4, cb5, cb6))
      cb.ax.set_ylabel(ylabels[i]; fontsize)
      cb.outline.set_linewidth(1.0)
   end

   return fig, cs, range1, range2
end

@everywhere function update_plot!(fig, cs, range1, range2, outdir, file)
   isfile(outdir*file[end-8:end-5]*".png") && return

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

   savefig(outdir*file[end-8:end-5]*".png", bbox_inches="tight")
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status, outdir, cmaps, norms, ticks, pArgs1, extent)
   fig, cs, range1, range2 = init_figure(cmaps, norms, ticks, pArgs1, extent)
   while true
      file = take!(jobs)
      update_plot!(fig, cs, range1, range2, outdir, file)
      put!(status, true)
   end
   close(fig)
end

################################################################################
files = glob("bulk*.vlsv", ".")

nfile = length(files)
# Set output directory
const outdir = "contour/"

# Set colormaps for continuous and divergent data
const cmaps = (matplotlib.cm.turbo, matplotlib.cm.RdBu_r)
const axisunit = EARTH
const extent = [0., 20., -20., 20.] # [RE], default: [-Inf32, Inf32, -Inf32, Inf32]

# Plotting range for each variable
const ρmin, ρmax   = 0.1, 15.0     # [amu/cc]
const vxmin, vxmax = -650.0, 650.0 # [km/s]
const vzmin, vzmax = -500.0, 500.0 # [km/s]
const pmin, pmax   = 0.0, 3.6      # [nPa]
const bmin, bmax   = -60.0, 60.    # [nT]
const emin, emax   = 0.0, 20.      # [mV/m]

meta = load(files[1])

pArgs1 = Vlasiator.set_args(meta, "proton/vg_rho", axisunit; normal=:none)

norm1, ticks1 = set_colorbar(Linear, ρmin, ρmax)
norm2, ticks2 = set_colorbar(Linear, vxmin, vxmax)
norm3, ticks3 = set_colorbar(Linear, vzmin, vzmax)
norm4, ticks4 = set_colorbar(Linear, pmin, pmax)
norm5, ticks5 = set_colorbar(Linear, bmin, bmax)
norm6, ticks6 = set_colorbar(Linear, emin, emax)

const norms = (norm1, norm2, norm3, norm4, norm5, norm6)
const ticks = (ticks1, ticks2, ticks3, ticks4, ticks5, ticks6)

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

@broadcast begin
   const fontsize = "x-large"
end

println("Total number of files: $nfile")
println("Running with $(nworkers()) workers...")

@async make_jobs(files) # Feed the jobs channel with all files to process.

@sync for p in workers()
   @async remote_do(do_work, p, jobs, status, outdir, cmaps, norms, ticks, pArgs1, extent)
end

let n = nfile
   t = @elapsed while n > 0 # wait for all jobs to complete
      take!(status)
      n -= 1
   end
   println("Finished in $(round(t, digits=2))s.")
end
```