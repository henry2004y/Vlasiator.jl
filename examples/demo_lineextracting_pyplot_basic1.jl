# Sample postprocessing script for extracting data alone a line across frames.
# Note: this script can be faster by replacing data instead of replotting everytime!
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob, PyPlot, Printf
using Vlasiator: RE # Earth radius, [m]

function main()
   files = glob("run4/bulk*.vlsv")

   point1 = [0e8, 0, 0]
   point2 = [1.9e8, 0, 0]

   cellids, distances, coords =
      let meta = load(files[1])
         getcellinline(meta, point1, point2)
      end

   # time density temperature vx
   inputs = [
   0.0 2.0e6 1.0e5 -5.0e5;
   99.0 2.0e6 1.0e5 -5.0e5;
   100.0 2.0e6 1.0e5 -8.0e5;
   200.0 2.0e6 1.0e5 -8.0e5;
   201.0 1.0e6 1.0e5 -8.0e5;
   499.0 1.0e6 1.0e5 -8.0e5;
   500.0 2.0e6 1.0e5 -5.0e5]

   ndigits = 4

   lim_rho = [0.0, 4.0e6]
   lim_v   = [-2e5, 10e5]

   fig, ax = plt.subplots(2,1, figsize=(12,5))

   for (i, file) in enumerate(files)
      fileout = lpad(i, ndigits, '0')*".png"

      meta = load(file)
      rho = readvariable(meta, "proton/vg_rho", cellids)
      v = readvariable(meta, "proton/vg_v", cellids)
      vx = Vector{Float64}(undef, size(v,1))

      for k in axes(v,1)
         vx[k] = v[k][1]
      end

      t = readparameter(meta, "time")

      ax[1].plot(coords[1,:] ./ RE, rho, label="rho")
      ax[1].set_ylim(lim_rho)
      ax[1].legend()
      ax[1].minorticks_on()

      ax[2].plot(coords[1,:] ./ RE, abs.(vx), label="vx")
      ax[2].set_ylim(lim_v)
      ax[2].legend()
      ax[2].minorticks_on()

      xlabel("X")

      j = searchsortedfirst(inputs[:,1], t) - 1
      if j == 0; j = 1; end
      str_title = 
         @sprintf "t=%4.1fs, rho=%4.1f/cc, vx=%4.1fkm/s" t inputs[j,2]/1e6 inputs[j,4]/1e3
      fig.suptitle(str_title, fontsize=14)

      savefig(fileout, bbox_inches="tight")

      ax[1].cla()
      ax[2].cla()
   end

   close(fig)
end

main()