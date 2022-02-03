# Sample script for plotting fieldlines with handpicked seeds, multi-processing version.
#
# To run on a single node,
# julia -p $ncores demo_fieldline_mp_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, PyCall, Printf, LaTeXStrings, FieldTracer
@everywhere using Vlasiator: set_args, prep2d, set_colorbar, RE

@everywhere function init_figure()
   fig, ax = plt.subplots(1, 1; num=myid(),
      figsize=(9, 9), constrained_layout=true)

   fontsize = "x-large"

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

   ax.set_xlabel(L"X [$R_E$]"; fontsize)
   ax.set_ylabel(L"Y [$R_E$]"; fontsize)

   ax.set_title("Density"; fontsize)

   x, y = Vlasiator.get_axis(pArgs)

   fakedata = zeros(Float32, length(y), length(x))
   c = ax.pcolormesh(x, y, fakedata; norm, cmap=matplotlib.cm.turbo)

   format = matplotlib.ticker.FormatStrFormatter("%.1f")
   cb1 = colorbar(c; ax, ticks, format)
   cb1.ax.set_ylabel("[amu/cc]"; fontsize)
   cb1.outline.set_linewidth(1.0)

   nseeds = 10

   fakeline = [0.0, 1.0]
   ls = [ax.plot(fakeline, fakeline,  color="w") for _ in 1:nseeds]

   return fig, ax, c, ls
end

@everywhere function process(fig, ax, c, ls, file)
   isfile("out/"*file[end-8:end-5]*".png") && return

   println("file = $file")
   meta = load(file)

   data = prep2d(meta, "proton/vg_rho", :mag)'
   c.set_array(data ./ 1f6)

   str_title = @sprintf "Density pulse run, t= %4.1fs" meta.time
   ax.set_title(str_title; fontsize)

   (;coordmin, coordmax, ncells) = meta
   dim_ = (1,2)

   # regular Cartesian mesh
   grid1 = range(coordmin[dim_[1]], coordmax[dim_[1]], length=ncells[dim_[1]]) 
   grid2 = range(coordmin[dim_[2]], coordmax[dim_[2]], length=ncells[dim_[2]])

   nseeds = length(ls)
   seeds = Matrix{Float64}(undef, 2, nseeds)
   for i in 1:nseeds
      seeds[1,i] = coordmin[dim_[1]] +
         (coordmax[dim_[1]] - coordmin[dim_[1]]) / nseeds * (i - 1)
      seeds[2,i] = -20RE
   end

   b = meta["vg_b_vol"]
   b1 = reshape(b[dim_[1],:], ncells[dim_[1]], ncells[dim_[2]])
   b2 = reshape(b[dim_[2],:], ncells[dim_[1]], ncells[dim_[2]])
   
   annotations = [child for child in ax.get_children() if
      pybuiltin(:isinstance)(child, matplotlib.text.Annotation)]
   
   for a in annotations
      a.remove()
   end

   for i = 1:nseeds
      startx, starty = seeds[:,i]
      x1, y1 = trace(b1, b2, startx, starty, grid1, grid2;
         ds=0.5, maxstep=4000, gridtype="ndgrid")
      x1 ./= RE
      y1 ./= RE
      if length(x1) < 5; continue; end
      ls[i][1].set_xdata(x1)
      ls[i][1].set_ydata(y1)
      add_arrow(ls[i][1])
   end

   savefig("out/"*file[end-8:end-5]*".png", bbox_inches="tight")
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status)
   fig, ax, c, ls = init_figure()
   while true
      file = take!(jobs)
      process(fig, ax, c, ls, file)
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
   colorscale = Linear
   axisunit = RE

   # Upper/lower limits for each variable
   const ρmin, ρmax = 0.0, 11.0 # [amu/cc]

   meta = load(files[1])

   const pArgs = set_args(meta, "proton/vg_rho", axisunit; normal=:none)
   const norm, ticks = set_colorbar(colorscale, ρmin, ρmax)
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