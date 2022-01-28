# Combined 1D/2D plots across multiple frames.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Glob, Printf, LaTeXStrings
using Vlasiator: set_args, prep2d, set_colorbar, set_plot

@assert VERSION ≥ v"1.7.0" "Compatible with Julia v1.7+!"
@assert matplotlib.__version__ >= "3.4" "Require Matplotlib version 3.4+ to use subfigure!"


struct Varminmax{T}
   "Density, [amu/cc]"
   ρmin::T
   ρmax::T
   "Velocity, [km/s]"
   vmin::T
   vmax::T
   "Pressure, [nPa]"
   pmin::T
   pmax::T
   "Magnetic field, [nT]"
   bmin::T
   bmax::T
   "Electric field, [nT]"
   emin::T
   emax::T
   "Alfven speed, [km/s]"
   vamin::T
   vamax::T
   "Sonic speed, [km/s]"
   vsmin::T
   vsmax::T
end

function init_figure(varminmax, loc, pArgs)
   (; ρmin, ρmax, vmin, vmax, pmin, pmax, bmin, bmax, emin, emax,
      vamin, vamax, vsmin, vsmax) = varminmax

   fig = plt.figure(constrained_layout=true, figsize=(12, 7.2))
   subfigs = fig.subfigures(1, 2, wspace=0.01, width_ratios=[2,1])

   axsL = subfigs[1].subplots(5, 1, sharex=true)
   axsR = subfigs[2].subplots(2, 1, sharex=true)

   # Set line plots' axes
   axsL[end].set_xlim(loc[1], loc[end])

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

   x, y = Vlasiator.get_axis(pArgs)
   fakedata = fill(NaN32, length(y), length(x))

   cnorm1, cticks1 = set_colorbar(Linear, vamin, vamax)
   cnorm2, cticks2 = set_colorbar(Linear, vsmin, vsmax)
   cmap = matplotlib.cm.turbo

   c1 = axsR[1].pcolormesh(x, y, fakedata; norm=cnorm1, cmap)
   c2 = axsR[2].pcolormesh(x, y, fakedata; norm=cnorm2, cmap)

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


"Update frame."
function process(plotargs, file)
   isfile("../out/"*file[end-8:end-5]*".png") && return

   fig, subfigs, ls, vlines, cs = plotargs

   meta = load(file)

   p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(x -> x*x, v_extract, dims=1) |> vec
   pram_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]

   bz = readvariable(meta, "vg_b_vol", cellids)[3,:] .* 1e9 #[nT]

   ls[1][1].set_ydata(rho_extract ./ 1e6)
   ls[2][1].set_ydata(v_extract[1,:] ./ 1e3)
   ls[3][1].set_ydata(pram_extract)
   ls[4][1].set_ydata(p_extract)
   ls[5][1].set_ydata(bz)

   imagnetopause_ = findfirst(<(0.0), bz)
   for vline in vlines
      update_vline(vline, loc[imagnetopause_])
   end

   str_title = @sprintf "Sun-Earth line, t= %4.1fs" meta.time
   subfigs[1].suptitle(str_title, fontsize=14)

   data = prep2d(meta, "VA", :z)'
   cs[1].update(Dict("array" => data ./ 1e3))

   data = prep2d(meta, "VS", :z)'
   cs[2].update(Dict("array" => data ./ 1e3))

   savefig("../out/"*file[end-8:end-5]*".png", bbox_inches="tight")
   return
end

function update_vline(h, x)
   seg_old = h.get_segments()
   ymin = seg_old[1][1, 2]
   ymax = seg_old[1][2, 2]

   seg_new = [[x ymin; x ymax]]

   h.set_segments(seg_new)
end

####### Main

files = glob("bulk*.vlsv", "../run_rho2_bz-5_timevarying_startfrom300s")
nfile = length(files)

meta = load(files[1])

x1, x2 = 7.0, 18.0 # Earth radii
point1 = [x1, 0, 0] .* Vlasiator.Re
point2 = [x2, 0, 0] .* Vlasiator.Re

cellids, distances, coords = getcellinline(meta, point1, point2)

loc = range(x1, x2, length=length(cellids))

pArgs = set_args(meta, "fakename", RE; normal=:none)

# Upper/lower limits for each variable
ρmin, ρmax = 0.0, 10.0     # [amu/cc]
vmin, vmax = -640.0, 100.0 # [km/s]
pmin, pmax = 0.0, 1.82     # [nPa]
bmin, bmax = -25.0, 60.0   # [nT]
emin, emax = -5.0, 5.0     # [mV/m]
vamin, vamax = 0.0, 250.0  # [km/s]
vsmin, vsmax = 0.0, 400.0  # [km/s]

varminmax =
   Varminmax(ρmin, ρmax, vmin, vmax, pmin, pmax, bmin, bmax, emin, emax,
   vamin, vamax, vsmin, vsmax)


plotargs = init_figure(varminmax, loc, pArgs)

# Loop over snapshots
for (i, file) in enumerate(files)
   println("i = $i/$nfile, file = $(basename(file))")
   process(plotargs, file)
end

close(plotargs[1])
println("Finished!")