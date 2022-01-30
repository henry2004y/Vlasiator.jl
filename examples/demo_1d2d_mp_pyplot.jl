# Combined 1D/2D plots across multiple frames, multi-process version.
#
# To run on a single node,
# julia -p $ncores demo_1d2d_mp_pyplot.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Distributed, ParallelDataTransfer, Glob
@everywhere using Vlasiator, PyPlot, Printf, LaTeXStrings
@everywhere using Vlasiator: set_args, prep2d, set_colorbar

@assert matplotlib.__version__ ≥ "3.4" "Require Matplotlib version 3.4+ to use subfigure!"

@everywhere function init_figure(x1, x2)
   fig = plt.figure(myid(), constrained_layout=true, figsize=(12, 7.2))
   subfigs = fig.subfigures(1, 2, wspace=0.01, width_ratios=[2,1])

   axsL = subfigs[1].subplots(5, 1, sharex=true)
   axsR = subfigs[2].subplots(2, 1, sharex=true)

   # Set line plots' axes
   axsL[end].set_xlim(x1, x2)

   axsL[1].set_ylim(ρmin, ρmax)
   axsL[2].set_ylim(vmin, vmax)
   axsL[3].set_ylim(pmin, pmax)
   axsL[4].set_ylim(bmin, bmax)
   axsL[5].set_ylim(emin, emax)
   for ax in axsL
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      ax.grid(true)
   end

   axsL[end].set_xlabel(L"x [$R_E$]"; fontsize=14)
   axsL[1].set_ylabel("n [amu/cc]";   fontsize=14)
   axsL[2].set_ylabel("V [km/s]";     fontsize=14)
   axsL[3].set_ylabel("P [nPa]";      fontsize=14)
   axsL[4].set_ylabel("B [nT]";       fontsize=14)
   axsL[5].set_ylabel("E [mV/m]";     fontsize=14)

   fakeline = loc
   l1 = axsL[1].plot(loc, fakeline, label="Proton density", color="#1f77b4")
   l2 = axsL[2].plot(loc, fakeline, label="Vx",             color="#1f77b4")
   l3 = axsL[2].plot(loc, fakeline, label="Vy",             color="#ff7f0e")
   l4 = axsL[2].plot(loc, fakeline, label="Vz",             color="#2ca02c")
   l5 = axsL[3].plot(loc, fakeline, label="Ram",            color="#1f77b4")
   l6 = axsL[3].plot(loc, fakeline, label="Thermal",        color="#ff7f0e")
   l7 = axsL[4].plot(loc, fakeline, label="Bx",             color="#1f77b4")
   l8 = axsL[4].plot(loc, fakeline, label="By",             color="#ff7f0e")
   l9 = axsL[4].plot(loc, fakeline, label="Bz",             color="#2ca02c")
   l10= axsL[5].plot(loc, fakeline, label="Ex",             color="#1f77b4")
   l11= axsL[5].plot(loc, fakeline, label="Ey",             color="#ff7f0e")
   l12= axsL[5].plot(loc, fakeline, label="Ez",             color="#2ca02c")

   ls = (l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12)

   axsL[2].legend(;loc="lower left",  ncol=3, frameon=false, fontsize=12)
   axsL[3].legend(;loc="upper left",  ncol=2, frameon=false, fontsize=12)
   axsL[4].legend(;loc="upper right", ncol=3, frameon=false, fontsize=12)
   axsL[5].legend(;loc="lower right", ncol=3, frameon=false, fontsize=12)

   vl1 = axsL[1].vlines(loc[1], ρmin, ρmax; colors="r", linestyle="dashed", alpha=0.5)
   vl2 = axsL[2].vlines(loc[1], vmin, vmax; colors="r", linestyle="dashed", alpha=0.5)
   vl3 = axsL[3].vlines(loc[1], pmin, pmax; colors="r", linestyle="dashed", alpha=0.5)
   vl4 = axsL[4].vlines(loc[1], bmin, bmax; colors="r", linestyle="dashed", alpha=0.5)
   vl5 = axsL[5].vlines(loc[1], emin, emax; colors="r", linestyle="dashed", alpha=0.5)

   hl2 = axsL[2].hlines(0.0, loc[1], loc[end]; colors="k", linestyle="dashed", alpha=0.2)
   hl4 = axsL[4].hlines(0.0, loc[1], loc[end]; colors="k", linestyle="dashed", alpha=0.2)
   hl5 = axsL[5].hlines(0.0, loc[1], loc[end]; colors="k", linestyle="dashed", alpha=0.2)

   vlines = (vl1, vl2, vl3, vl4, vl5)

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

      ax.set_ylabel(L"Y [$R_E$]"; fontsize=14)
   end

   axsR[2].set_xlabel(L"X [$R_E$]"; fontsize=14)

   axsR[1].set_title("Alfven speed", fontsize=14)
   axsR[2].set_title("Sound speed", fontsize=14)

   x, y = Vlasiator.get_axis(pArgs1)
   fakedata = fill(NaN32, length(y), length(x))

   c1 = axsR[1].pcolormesh(x, y, fakedata, norm=cnorm1, cmap=cmap)
   c2 = axsR[2].pcolormesh(x, y, fakedata, norm=cnorm2, cmap=cmap)

   rInner = 31.8e6 # [m]
   circle1 = plt.Circle((0, 0), rInner/Vlasiator.Re, color="w")
   circle2 = plt.Circle((0, 0), rInner/Vlasiator.Re, color="w")
   axsR[1].add_patch(circle1)
   axsR[2].add_patch(circle2)

   cb1 = colorbar(c1; ax=axsR[1], ticks=cticks1, fraction=0.046, pad=0.04, extend="max")
   cb1.ax.set_ylabel("[km/s]"; fontsize=14)

   cb2 = colorbar(c2; ax=axsR[2], ticks=cticks2, fraction=0.046, pad=0.04)
   cb2.ax.set_ylabel("[km/s]"; fontsize=14)

   #fig.suptitle("Density Pulse Run", fontsize="xx-large")

   cs = (c1, c2)

   return fig, subfigs, ls, vlines, cs
end

@everywhere function update_vline(h, x)
   seg_old = h.get_segments()
   ymin = seg_old[1][1, 2]
   ymax = seg_old[1][2, 2]

   seg_new = [[x ymin; x ymax]]

   h.set_segments(seg_new)
end

@everywhere function process(subfigs, ls, vlines, cs, file, cellids)
   isfile("../out/"*file[end-8:end-5]*".png") && return

   println("file = $(basename(file))")
   meta = load(file)

   rho = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v = readvariable(meta, "proton/vg_v", cellids)
   p = readvariable(meta, "vg_pressure", cellids) .* 1f9 |> vec # [nPa]

   vmag2 = sum(x -> x*x, v, dims=1) |> vec
   pram = rho .* Vlasiator.mᵢ .* vmag2 .* 1f9 # [nPa]

   b = readvariable(meta, "vg_b_vol", cellids) .* 1f9 #[nT]
   e = readvariable(meta, "vg_e_vol", cellids) .* 1f3 #[mV/m]

   ls[1][1].set_ydata(rho ./ 1f6)
   ls[2][1].set_ydata(@views v[1,:] ./ 1f3)
   ls[3][1].set_ydata(@views v[2,:] ./ 1f3)
   ls[4][1].set_ydata(@views v[3,:] ./ 1f3)
   ls[5][1].set_ydata(pram)
   ls[6][1].set_ydata(p)
   ls[7][1].set_ydata(@view b[1,:])
   ls[8][1].set_ydata(@view b[2,:])
   ls[9][1].set_ydata(@view b[3,:])
   ls[10][1].set_ydata(@view e[1,:])
   ls[11][1].set_ydata(@view e[2,:])
   ls[12][1].set_ydata(@view e[3,:])

   imagnetopause_ = findfirst(<(0.0), @views b[3,:])
   for vline in vlines
      update_vline(vline, loc[imagnetopause_])
   end

   str_title = @sprintf "Sun-Earth line, t= %4.1fs" meta.time
   subfigs[1].suptitle(str_title, fontsize=14)

   data = prep2d(meta, "VA", :z)'
   cs[1].set_array(data ./ 1f3)

   data = prep2d(meta, "VS", :z)'
   cs[2].set_array(data ./ 1f3)

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
files = glob("bulk*.vlsv", ".")
const nfile = length(files)
@passobj 1 workers() files

@broadcast begin
   # Set contour plots' axes and colorbars
   const cmap = matplotlib.cm.turbo
   axisunit = RE
   # Upper/lower limits for each variable
   const ρmin, ρmax = 0.0, 10.0     # [amu/cc]
   const vmin, vmax = -640.0, 100.0 # [km/s]
   const pmin, pmax = 0.0, 1.82     # [nPa]
   const bmin, bmax = -25.0, 60.0   # [nT]
   const emin, emax = -5.0, 5.0     # [mV/m]
   const vamin, vamax = 50.0, 350.0 # [km/s]
   const vsmin, vsmax = 50.0, 350.0 # [km/s]

   meta = load(files[1])

   const pArgs1 = set_args(meta, "VA", axisunit; normal=:none)
   const cnorm1, cticks1 = set_colorbar(Linear, vamin, vamax)
   const cnorm2, cticks2 = set_colorbar(Linear, vsmin, vsmax)
end

const jobs   = RemoteChannel(()->Channel{String}(nfile))
const status = RemoteChannel(()->Channel{Bool}(nworkers()))

x1, x2 = 7.0, 20.0 # Earth radii
point1 = [x1, 0, 0] .* Vlasiator.Re
point2 = [x2, 0, 0] .* Vlasiator.Re

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