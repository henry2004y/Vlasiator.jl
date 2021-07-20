# Pressure plot alone a line.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

filenames = glob("bulk*.vlsv")
#filenames = ["bulk.0001020.vlsv"]

meta = load(filenames[1])

Re = Vlasiator.Re # Earth radii
x1, x2 = 5.0, 29.0
point1 = [x1, 0, 0] .* Re
point2 = [x2, 0, 0] .* Re

cellids, distances, coords = getcellinline(meta, point1, point2)

close(meta.fid)

fig, ax = plt.subplots()

for fname in filenames
   local meta = load(fname)

   p_extract = readvariable(meta, "vg_pressure", cellids) |> vec
   rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
   v_extract = readvariable(meta, "proton/vg_v", cellids)
   vmag2_extract = sum(v_extract.^2, dims=1) |> vec
   pdyn_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract
   loc = range(x1, x2, length=length(rho_extract))

   ax.plot(loc, pdyn_extract, label="dynamic")
   ax.plot(loc, p_extract, label="thermal")

   xlim(x1, x2)
   ylim(0.0, 1.9e-9)
   grid()

   xlabel("x [Re]")
   ylabel("Pressure [Pa]")
   title("t = $(meta.time)")
   legend()
   savefig("p_yz0line_510s.png",bbox_inches="tight")
   ax.cla()
end