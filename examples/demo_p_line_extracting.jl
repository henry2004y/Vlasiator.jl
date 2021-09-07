# Pressure plot alone a line across multiple frames.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Glob, Printf

filenames = glob("bulk*.vlsv", "run_rho2_bz-5_timevarying_startfrom300s")
nfile = length(filenames)

meta = load(filenames[1])

Re = Vlasiator.Re # Earth radii
x1, x2 = 8.0, 29.0
point1 = [x1, 0, 0] .* Re
point2 = [x2, 0, 0] .* Re

cellids, distances, coords = getcellinline(meta, point1, point2)

close(meta.fid)

fig, ax = plt.subplots(figsize=(8, 4.8))

for (i, fname) in enumerate(filenames)
   println("i = $i/$nfile, filename = $fname")
   local meta = load(fname)

   p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(v_extract.^2, dims=1) |> vec
   pdyn_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]
   loc = range(x1, x2, length=length(rho_extract))

   ax.plot(loc, pdyn_extract, label="dynamic")
   ax.plot(loc, p_extract, label="thermal")

   ax.set_xlim(x1, x2)
   ax.set_ylim(0.0, 1.8)
   ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
   grid()

   xlabel("x [Re]", fontsize=14)
   ylabel("Pressure [nPa]", fontsize=14)
   str_title = @sprintf "t= %4.1fs" meta.time
   title(str_title, fontsize=14)
   legend(loc="upper left", fontsize=14)
   savefig("out/"*lpad(i, 4, '0')*".png", bbox_inches="tight")
   ax.cla()
end

close(fig)
println("plot finished!")