# Combined 1D/2D plots across multiple frames, multi-process version.
#
# To run on a single node,
# julia -p $ncores demo_1d2d_parallel_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, Printf, LaTeXStrings
@everywhere using Vlasiator: set_args, plot_prep2d, set_colorbar, set_plot

@assert matplotlib.__version__ ≥ "3.4" "Require Matplotlib version 3.4+ to use subfigure!"

@everywhere function init_figure(x1, x2)
   fig = plt.figure(myid(), constrained_layout=true, figsize=(12, 12))
   subfigs = fig.subfigures(1, 2, wspace=0.05)

   axsL = subfigs[1].subplots(4, 1, sharex=true)
   axsR = subfigs[2].subplots(2, 1, sharex=true, sharey=true)

   # Set line plots' axes
   axsL[end].set_xlim(x1, x2)

   axsL[1].set_ylim(ρmin, ρmax)
   axsL[2].set_ylim(vmin, vmax)
   axsL[3].set_ylim(pmin, pmax)
   axsL[4].set_ylim(bmin, bmax)
   for ax in axsL
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.grid(true)
   end

   axsL[end].set_xlabel(L"x [$R_E$]"; fontsize)
   axsL[1].set_ylabel("Density [amu/cc]";    fontsize)
   axsL[2].set_ylabel("Velocity [km/s]";     fontsize)
   axsL[3].set_ylabel("Pressure [nPa]";      fontsize)
   axsL[4].set_ylabel("Magnetic field [nT]"; fontsize)

   fakeline = loc
   l1 = axsL[1].plot(loc, fakeline, label="Proton density", color="#1f77b4")
   l2 = axsL[2].plot(loc, fakeline, label="Vx",             color="#1f77b4")
   l3 = axsL[3].plot(loc, fakeline, label="Dynamic",        color="#1f77b4")
   l4 = axsL[3].plot(loc, fakeline, label="Thermal",        color="#ff7f0e")
   l5 = axsL[4].plot(loc, fakeline, label="Bz",             color="#1f77b4")

   ls = (l1, l2, l3, l4, l5)

   axsL[2].legend(;loc="upper right", fontsize)
   axsL[3].legend(;loc="lower right", fontsize)
   axsL[4].legend(;loc="upper right", fontsize)

   vl1 = axsL[1].vlines(loc[imagnetopause_], ρmin, ρmax;
      colors="r", linestyle="dashed", alpha=0.5)

   vl2 = axsL[2].vlines(loc[imagnetopause_], vmin, vmax;
      colors="r", linestyle="dashed", alpha=0.5)

   vl3 = axsL[3].vlines(loc[imagnetopause_], pmin, pmax;
      colors="r", linestyle="dashed", alpha=0.5)

   vl4 = axsL[4].vlines(loc[imagnetopause_], bmin, bmax;
       colors="r", linestyle="dashed", alpha=0.5)
   hl4 = axsL[4].hlines(0.0, loc[1], loc[end]; colors="k", linestyle="dashed", alpha=0.2)

   vlines = (vl1, vl2, vl3, vl4)

   for ax in axsR
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

      ax.set_ylabel(L"Y [$R_E$]"; fontsize)
   end

   axsR[2].set_xlabel(L"X [$R_E$]"; fontsize)

   axsR[1].set_title("Alfven speed", fontsize="x-large")
   axsR[2].set_title("Sound speed", fontsize="x-large")

   plotrange, sizes, axisunit = pArgs1.plotrange, pArgs1.sizes, pArgs1.axisunit
   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   fakedata = zeros(Float32, length(y), length(x))

   c1 = axsR[1].pcolormesh(x, y, fakedata, norm=cnorm1, cmap=cmap, shading="nearest")
   c2 = axsR[2].pcolormesh(x, y, fakedata, norm=cnorm2, cmap=cmap, shading="nearest")

   cb1 = colorbar(c1; ax=axsR[1], ticks=cticks1, fraction=0.046, pad=0.04)
   cb1.ax.set_ylabel("[km/s]"; fontsize)
   cb1.outline.set_linewidth(1.0)

   cb2 = colorbar(c2; ax=axsR[2], ticks=cticks2, fraction=0.046, pad=0.04)
   cb2.ax.set_ylabel("[km/s]"; fontsize)
   cb2.outline.set_linewidth(1.0)

   fig.suptitle("Density Pulse Run", fontsize="xx-large")

   cs = (c1, c2)

   return fig, subfigs, ls, vlines, cs
end

function update_vline(h, x)
   seg_old = h.get_segments()
   ymin = seg_old[1][1, 2]
   ymax = seg_old[1][2, 2]

   seg_new = [[[x, ymin] [x, ymax]]]

   h.set_segments(seg_new)
end

@everywhere function process(subfigs, ls, vlines, cs, file, cellids)
   isfile("../out/"*file[end-8:end-5]*".png") && return

   println("file = $(basename(file))")
   meta = load(file)

   p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(v_extract.^2, dims=1) |> vec
   pdyn_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]

   bz = readvariable(meta, "vg_b_vol", cellids)[3,:] .* 1e9 #[nT]

   ls[1][1].set_ydata(rho_extract ./ 1e6)
   ls[2][1].set_ydata(v_extract[1,:] ./ 1e3)
   ls[3][1].set_ydata(pdyn_extract)
   ls[4][1].set_ydata(p_extract)
   ls[5][1].set_ydata(bz)

   imagnetopause_ = findfirst(<(0.0), bz)
   for vline in vlines
      update_vline(vline, loc[imagnetopause_])
   end

   str_title = @sprintf "Sun-Earth line, t= %4.1fs" meta.time
   subfigs[1].suptitle(str_title, fontsize="x-large")

   _, _, data = plot_prep2d(meta, "VA", pArgs1, :z)
   cs[1].set_array(data ./ 1e3)

   _, _, data = plot_prep2d(meta, "VS", pArgs2, :z)
   cs[2].set_array(data ./ 1e3)

   savefig("../out/"*file[end-8:end-5]*".png", bbox_inches="tight")
   return
end

function make_jobs(files)
   for f in files
       put!(jobs, f)
   end
end

@everywhere function do_work(jobs, status)
   fig, subfigs, ls, vlines, cs = init_figure(x1, x2)
   while true
      file = take!(jobs)
      process(subfigs, ls, vlines, cs, file, cellids)
      put!(status, true)
   end
   close(fig)
end

############################################################################################
files = glob("bulk*.vlsv", "../run_rho2_bz-5_timevarying_startfrom300s")
const nfile = length(files)
@passobj 1 workers() files

@broadcast begin
   # Set contour plots' axes and colorbars
   const cmap = matplotlib.cm.turbo
   colorscale = Linear
   axisunit = RE
   # Upper/lower limits for each variable
   const ρmin, ρmax = 0.0, 10.0     # [amu/cc]
   const vmin, vmax = -640.0, 0.0   # [km/s]
   const pmin, pmax = 0.0, 1.82     # [nPa]
   const bmin, bmax = -25.0, 60.0   # [nT]
   const vamin, vamax = 0.0, 250.0  # [km/s]
   const vsmin, vsmax = 30.0, 350.0 # [km/s]

   meta = load(files[1])

   const pArgs1 = set_args(meta, "VA", axisunit, colorscale;
      normal=:none, vmin=vamin, vmax=vamax)
   const cnorm1, cticks1 = set_colorbar(pArgs1)

   const pArgs2 = set_args(meta, "VS", axisunit, colorscale;
      normal=:none, vmin=vsmin, vmax=vsmax)
   const cnorm2, cticks2 = set_colorbar(pArgs2)

   const fontsize = 14
end

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

Re = Vlasiator.Re # Earth radii
x1, x2 = 8.0, 29.0
point1 = [x1, 0, 0] .* Re
point2 = [x2, 0, 0] .* Re

meta = load(files[1])
cellids, _, _ = getcellinline(meta, point1, point2)

passobj(1, workers(), [:x1, :x2, :cellids])
@broadcast const loc = range(x1, x2, length=length(cellids))

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