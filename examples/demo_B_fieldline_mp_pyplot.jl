# Sample script for plotting fieldlines with handpicked seeds, multi-processing version.
#
# To run on a single node,
# julia -p $ncores demo_fieldline_mp_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, PyCall, Printf, LaTeXStrings, FieldTracer
@everywhere using Vlasiator: RE

function generate_seeds(coordmin, coordmax, dim_, nseeds)
   seeds = Matrix{Float64}(undef, 2, nseeds)
   for i in 1:nseeds
      seeds[1,i] = coordmin[dim_[1]] +
         (coordmax[dim_[1]] - coordmin[dim_[1]]) / nseeds * (i - 1)
      seeds[2,i] = -20RE
   end
   seeds
end

@everywhere function init_figure(pArgs, norm, ticks, seeds, extent)
   fig, ax = plt.subplots(1, 1; num=myid(),
      figsize=(6, 8), constrained_layout=true)

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

   ax.set_xlabel(pArgs.strx; fontsize)
   ax.set_ylabel(pArgs.stry; fontsize)

   ax.set_title("Density"; fontsize)

   x1, x2 = Vlasiator.get_axis(pArgs)

   range1 = searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2])
   range2 = searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])

   fakedata = zeros(Float32, length(range2), length(range1))
   c = ax.pcolormesh(x1[range1], x2[range2], fakedata; norm, cmap=matplotlib.cm.turbo)

   format = matplotlib.ticker.FormatStrFormatter("%.1f")
   cb1 = colorbar(c; ax, ticks, format)
   cb1.ax.set_ylabel("[amu/cc]"; fontsize)
   cb1.outline.set_linewidth(1.0)

   fakeline = [0.0, 1.0]
   ls = [ax.plot(fakeline, fakeline,  color="w") for _ in 1:size(seeds,2)]

   return fig, ax, c, ls, range1, range2
end

@everywhere function update_plot!(ax, c, ls, range1, range2, dim_, seeds, grid1, grid2,
   outdir, file)
   isfile(outdir*file[end-8:end-5]*".png") && return

   println("file = $file")
   meta = load(file)

   data = Vlasiator.prep2d(meta, "proton/vg_rho", :mag)'
   c.set_array(data[range2,range1] ./ 1f6)

   str_title = @sprintf "Density pulse run, t= %4.1fs" meta.time
   ax.set_title(str_title; fontsize="x-large")

   b = meta["vg_b_vol"]
   b1 = reshape(b[dim_[1],:], meta.ncells[dim_[1]], meta.ncells[dim_[2]])
   b2 = reshape(b[dim_[2],:], meta.ncells[dim_[1]], meta.ncells[dim_[2]])
   # Find existing arrow annotations
   annotations = [child for child in ax.get_children() if
      pybuiltin(:isinstance)(child, matplotlib.text.Annotation)]
   # Remove existing arrows
   for a in annotations
      a.remove()
   end
   # Add new arrows along field lines
   for i in axes(seeds,2)
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

   savefig(outdir*file[end-8:end-5]*".png", bbox_inches="tight")
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status,
   outdir, pArgs, norm, ticks, grid1, grid2, dim_, seeds, extent)

   fig, ax, c, ls, range1, range2 = init_figure(pArgs, norm, ticks, seeds, extent)

   while true
      file = take!(jobs)
      update_plot!(ax, c, ls, range1, range2, dim_, seeds, grid1, grid2, outdir, file)
      put!(status, true)
   end
   close(fig)
end

############################################################################################
files = glob("bulk*.vlsv", ".")

nfile = length(files)
# Set output directory
outdir = "out/"

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

axisunit = EARTH # contour plot axes unit
extent = [0., 20., -20., 20.] # [RE], default full domain: [-Inf32, Inf32, -Inf32, Inf32]

# Upper/lower limits for each variable
ρmin, ρmax = 0.0, 11.0      # [amu/cc]

meta = load(files[1])
# Construct pieces for plotting
pArgs = Vlasiator.set_args(meta, "proton/vg_rho", axisunit; normal=:none)
norm, ticks = Vlasiator.set_colorbar(Linear, ρmin, ρmax)

# Mark spatial dimensions
dim_ = pArgs.stry[1] == 'Z' ? (1,3) : (1,2)

(;coordmin, coordmax, ncells) = meta

# Generate regular Cartesian range
grid1 = range(coordmin[dim_[1]], coordmax[dim_[1]], length=ncells[dim_[1]]) 
grid2 = range(coordmin[dim_[2]], coordmax[dim_[2]], length=ncells[dim_[2]])
# Generate seeds for in-plane field line tracing
nseeds = 10
seeds = generate_seeds(coordmin, coordmax, dim_, nseeds)


println("Total number of files: $nfile")
println("Running with $(nworkers()) workers...")

@async make_jobs(files) # Feed the jobs channel with all files to process.

@sync for p in workers()
   @async remote_do(do_work, p, jobs, status,
      outdir, pArgs, norm, ticks, grid1, grid2, dim_, seeds, extent)
end

let n = nfile
   t = @elapsed while n > 0 # wait for all jobs to complete
      take!(status)
      n -= 1
   end
   println("Finished in $(round(t, digits=2))s.")
end