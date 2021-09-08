# Combined 1D/2D plots across multiple frames.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Glob, Printf, LaTeXStrings
using Vlasiator: set_args, plot_prep2d, set_colorbar, set_plot

@assert matplotlib.__version__ >= "3.4" "Matplotlib version too old: subfigure is not supported!"

filenames = glob("bulk*.vlsv", ".")
nfile = length(filenames)

meta = load(filenames[1])

Re = Vlasiator.Re # Earth radii
x1, x2 = 8.0, 29.0
point1 = [x1, 0, 0] .* Re
point2 = [x2, 0, 0] .* Re

cellids, distances, coords = getcellinline(meta, point1, point2)

loc = range(x1, x2, length=length(cellids))

close(meta.fid)

# Upper/lower limits for each variable
ρmin, ρmax = 0.0, 10.0     # [amu/cc]
vmin, vmax = -640.0, 0.0   # [km/s]
pmin, pmax = 0.0, 1.82     # [nPa]
bmin, bmax = -25.0, 60.0   # [nT]
vamin, vamax = 0.0, 250.0  # [km/s]
vsmin, vsmax = 0.0, 400.0  # [km/s]

fig = plt.figure(constrained_layout=true, figsize=(12, 12))
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

axsL[end].set_xlabel(L"x [$R_E$]", fontsize=14)
axsL[1].set_ylabel("Density [amu/cc]", fontsize=14)
axsL[2].set_ylabel("Velocity [km/s]", fontsize=14)
axsL[3].set_ylabel("Pressure [nPa]", fontsize=14)
axsL[4].set_ylabel("Magnetic field [nT]", fontsize=14)

# Set contour plots' axes and colorbars
cmap = matplotlib.cm.turbo
colorscale = Linear
axisunit = RE
op = :z
      
pArgs1 = set_args(meta, "VA", axisunit, colorscale; normal=:none, vmin=vamin, vmax=vamax)
cnorm1, cticks1 = set_colorbar(pArgs1)

pArgs2 = set_args(meta, "VS", axisunit, colorscale; normal=:none, vmin=vsmin, vmax=vsmax)
cnorm2, cticks2 = set_colorbar(pArgs2)

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

   ax.set_ylabel(L"Y [$R_E$]", fontsize=14)
end

axsR[2].set_xlabel(L"X [$R_E$]", fontsize=14)

axsR[1].set_title("Alfven speed", fontsize="x-large")
axsR[2].set_title("Sound speed", fontsize="x-large")

fig.suptitle("Density Pulse Run", fontsize="xx-large")


# Loop over snapshots
for (i, fname) in enumerate(filenames)
   println("i = $i/$nfile, filename = $fname")
   local meta = load(fname)

   p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(v_extract.^2, dims=1) |> vec
   pdyn_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]

   bz = readvariable(meta, "vg_b_vol", cellids)[3,:] .* 1e9 #[nT]

   imagnetopause_ = findfirst(<(0.0), bz)

   axsL[1].plot(loc, rho_extract ./ 1e6, label="Proton density", color="#1f77b4")
   axsL[1].vlines(loc[imagnetopause_], ρmin, ρmax, colors="r", linestyle="dashed", alpha=0.5)

   axsL[2].plot(loc, v_extract[1,:] ./ 1e3, label="Vx", color="#1f77b4")
   axsL[2].vlines(loc[imagnetopause_], vmin, vmax, colors="r", linestyle="dashed", alpha=0.5)

   axsL[3].plot(loc, pdyn_extract, label="Dynamic", color="#1f77b4")
   axsL[3].plot(loc, p_extract, label="Thermal", color="#ff7f0e")
   axsL[3].vlines(loc[imagnetopause_], pmin, pmax, colors="r", linestyle="dashed", alpha=0.5)

   axsL[4].plot(loc, bz, label="Bz", color="#1f77b4")
   axsL[4].hlines(0.0, loc[1], loc[end], colors="k", linestyle="dashed", alpha=0.2)
   axsL[4].vlines(loc[imagnetopause_], bmin, bmax, colors="r", linestyle="dashed", alpha=0.5)

   str_title = @sprintf "Sun-Earth line, t= %4.1fs" meta.time
   subfigs[1].suptitle(str_title, fontsize="x-large")

   axsL[2].legend(loc="upper right", fontsize=14)
   axsL[3].legend(loc="lower right", fontsize=14)
   axsL[4].legend(loc="upper right", fontsize=14)

   x, y, data = plot_prep2d(meta, "VA", pArgs1, op, axisunit) 
   c1 = axsR[1].pcolormesh(x, y, data ./ 1e3, norm=cnorm1, cmap=cmap, shading="auto")

   x, y, data = plot_prep2d(meta, "VS", pArgs2, op, axisunit) 
   c2 = axsR[2].pcolormesh(x, y, data ./ 1e3, norm=cnorm2, cmap=cmap, shading="auto")

   if i == 1
      cb1 = colorbar(c1; ax=axsR[1], ticks=cticks1, fraction=0.046, pad=0.04)
      cb1_title = cb1.ax.set_ylabel("[km/s]", fontsize=14)
      cb1.outline.set_linewidth(1.0)

      cb2 = colorbar(c2; ax=axsR[2], ticks=cticks2, fraction=0.046, pad=0.04)
      cb2_title = cb2.ax.set_ylabel("[km/s]", fontsize=14)
      cb2.outline.set_linewidth(1.0)
   end

   savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")

   for ax in axsL
      for line in ax.get_lines()
         line.remove()
      end
   end
end

close(fig)
println("Finished!")


