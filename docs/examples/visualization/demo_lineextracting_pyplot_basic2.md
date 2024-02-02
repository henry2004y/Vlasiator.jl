# ---
# title: Plotting variable along a line
# id: demo_plot_line
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot variable along a line
# ---

This demo shows how to plot pressure along a line across multiple frames.
```julia
using Vlasiator, PyPlot, Printf
using Vlasiator: RE # Earth radius [m]

function main()
   files = filter(contains(r"^bulk.*\.vlsv$"), readdir())
   nfile = length(files)

   x1, x2 = 8.0, 29.0
   point1 = [x1, 0, 0] .* RE
   point2 = [x2, 0, 0] .* RE

   cellids, distances, coords =
      load(files[1]) do meta
         getcellinline(meta, point1, point2)
      end

   fig, ax = plt.subplots(figsize=(8, 4.8))

   for (i, file) in enumerate(files)
      println("i = $i/$nfile, file = $file")
      meta = load(file)

      p_extract = readvariable(meta, "vg_pressure", cellids) .* 1e9 |> vec # [nPa]
      rho_extract = readvariable(meta, "proton/vg_rho", cellids) |> vec
      v_extract = readvariable(meta, "proton/vg_v", cellids)
      vmag2_extract = sum(x -> x*x, v_extract, dims=1) |> vec
      pram_extract = rho_extract .* Vlasiator.mᵢ .* vmag2_extract .* 1e9 # [nPa]
      loc = range(x1, x2, length=length(rho_extract))

      ax.plot(loc, pram_extract, label="ram")
      ax.plot(loc, p_extract, label="thermal")

      ax.set_xlim(x1, x2)
      ax.set_ylim(0.0, 1.8)
      ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
      grid(true)

      xlabel("x [Re]", fontsize=14)
      ylabel("Pressure [nPa]", fontsize=14)
      str_title = @sprintf "t= %4.1fs" meta.time
      title(str_title, fontsize=14)
      legend(loc="upper left", fontsize=14)
      savefig("out/"*meta.name[end-8:end-5]*".png", bbox_inches="tight")
      ax.cla()
   end

   close(fig)
   println("plot finished!")
end

main()
```