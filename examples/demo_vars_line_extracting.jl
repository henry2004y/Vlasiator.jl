# Multiple variables alone a line across multiple frames.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Glob, Printf

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

ρmin, ρmax = 0.0, 10.0     # [amu/cc]
vmin, vmax = -640.0, 0.0   # [km/s]
pmin, pmax = 0.0, 1.82     # [nPa]
bmin, bmax = -25.0, 60.0   # [nT]

fontsize = 14

fig, axs = plt.subplots(4, 1, figsize=(10, 15), sharex=true, constrained_layout=true)

axs[end].set_xlim(x1, x2)

axs[1].set_ylim(ρmin, ρmax)
axs[2].set_ylim(vmin, vmax)
axs[3].set_ylim(pmin, pmax)
axs[4].set_ylim(bmin, bmax)
for ax in axs
   ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   ax.grid(true)
end

axs[end].set_xlabel("x [Re]"; fontsize)
axs[1].set_ylabel("Density [amu/cc]"; fontsize)
axs[2].set_ylabel("Velocity [km/s]"; fontsize)
axs[3].set_ylabel("Pressure [nPa]"; fontsize)
axs[4].set_ylabel("Magnetic field [nT]"; fontsize)

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

   axs[1].plot(loc, rho_extract ./ 1e6, label="Proton density", color="#1f77b4")
   vl1 = axs[1].vlines(loc[imagnetopause_], ρmin, ρmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[2].plot(loc, v_extract[1,:] ./ 1e3, label="Vx", color="#1f77b4")
   vl2 = axs[2].vlines(loc[imagnetopause_], vmin, vmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[3].plot(loc, pdyn_extract, label="Dynamic", color="#1f77b4")
   axs[3].plot(loc, p_extract, label="Thermal", color="#ff7f0e")
   vl3 = axs[3].vlines(loc[imagnetopause_], pmin, pmax;
      colors="r", linestyle="dashed", alpha=0.5)

   axs[4].plot(loc, bz, label="Bz", color="#1f77b4")
   hl4 = axs[4].hlines(0.0, loc[1], loc[end], colors="k", linestyle="dashed", alpha=0.2)
   vl4 = axs[4].vlines(loc[imagnetopause_], bmin, bmax;
       colors="r", linestyle="dashed", alpha=0.5)

   str_title = @sprintf "t= %4.1fs" meta.time
   axs[1].set_title(str_title; fontsize)

   axs[2].legend(;loc="upper right", fontsize)
   axs[3].legend(;loc="upper left",  fontsize)
   axs[4].legend(;loc="upper right", fontsize)

   savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")

   for ax in axs
      for line in ax.get_lines()
         line.remove()
      end
   end
   vl1.remove()
   vl2.remove()
   vl3.remove()
   vl4.remove()
   hl4.remove()
end

close(fig)
println("Finished!")