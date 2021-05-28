# Sample postprocessing script for extracting data alone a line across frames.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob, PyPlot, Printf

filenames = glob("run4/bulk*.vlsv")

meta = readmeta(filenames[1])

point1 = [0e8, 0, 0]
point2 = [1.9e8, 0, 0]

cellids, distances, coords = getcellinline(meta, point1, point2)

# time density temperature vx
inputs = [
0.0 2.0e6 1.0e5 -5.0e5;
99.0 2.0e6 1.0e5 -5.0e5;
100.0 2.0e6 1.0e5 -8.0e5;
200.0 2.0e6 1.0e5 -8.0e5;
201.0 1.0e6 1.0e5 -8.0e5;
499.0 1.0e6 1.0e5 -8.0e5;
500.0 2.0e6 1.0e5 -5.0e5]

nfiles = length(filenames)
ndigits = 4

lim_rho = [0.0, 4.0e6]
lim_v   = [-2e5, 10e5]

fig, ax = plt.subplots(2,1, figsize=(12,5))

for (i, filename) in enumerate(filenames)

   fnameout = lpad(i, ndigits, '0')*".png"

   local meta = readmeta(filename)
   local rho = readvariable(meta, "proton/vg_rho", cellids)
   local v = readvariable(meta, "proton/vg_v", cellids)
   local vx = Vector{Float64}(undef, size(v,1))

   for k = 1:size(v,1)
      vx[k] = v[k][1]
   end

   local t = readparameter(meta, "time")

   ax[1].plot(coords[1,:]./ Vlasiator.Re, rho, label="rho")
   ax[1].set_ylim(lim_rho)
   ax[1].legend()
   ax[1].minorticks_on()

   ax[2].plot(coords[1,:]./ Vlasiator.Re, abs.(vx), label="vx")
   ax[2].set_ylim(lim_v)
   ax[2].legend()
   ax[2].minorticks_on()

   xlabel("X")

   local j = searchsortedfirst(inputs[:,1], t) - 1
   if j == 0; j = 1; end
   local str_title = 
      @sprintf "t=%4.1fs, rho=%4.1f/cc, vx=%4.1fkm/s" t inputs[j,2]/1e6 inputs[j,4]/1e3
   fig.suptitle(str_title, fontsize=14)

   savefig(fnameout, bbox_inches="tight")

   ax[1].cla()
   ax[2].cla()
end

close()