# ---
# title: VDFs along a line
# id: demo_vdfs
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot velocity distribution functions
# ---

This demo shows how to plot velocity distribution functions along a line.

```julia
using VlasiatorPyPlot, JLD2
using Vlasiator: RE # Earth radius [m]

function main()
   ## Plot velocity space distributions
   xloc = range(7.0, 10.0, step=2.35)
   yloc = range(-0.0, 5.0, step=2.35)
   zloc = 0.0

   CIs = CartesianIndices((1:length(xloc), 1:length(yloc)))

   fig, axs = plt.subplots(length(xloc), length(yloc), sharex=true, sharey=true)

   for i in CIs
      loc = RE .* [xloc[i[1]], yloc[i[2]], zloc]
      vdfslice(meta, loc, axs[i]; verbose=false)
   end

   for a in axs[end,:] a.set_xlabel("vx [km/s]") end
   for a in axs[:  ,1] a.set_ylabel("vy [km/s]") end

   ## Subtract data of cells with VDF

   file = "test/data/bulk.1d.vlsv"
   meta = load(file)

   cells = getcellwithvdf(meta)

   ρ = zeros(Float64, size(cells))
   v = zeros(Float64, 3, size(cells)...)
   p = zeros(Float64, size(cells))

   vsize = meta.meshes["proton"].vblock_size .* meta.meshes["proton"].vblocks
   vcellids = [Int[] for _ in cells]
   vcellf = [Float64[] for _ in cells]

   for (i, cell) in enumerate(cells)
      ρ[i] = readvariable(meta, "proton/vg_rho", cell)[1]
      v[i] = readvariable(meta, "proton/vg_v", cell)[1]
      p[i] = readvariable(meta, "vg_pressure", cell)[1]

      vcellids[i], vcellf[i] = readvcells(meta, cell)
   end

   jldopen("example.jld2", "a+") do file
      file["file1/density"] = ρ
      file["file1/velocity"] = v
      file["file1/pressure"] = p
      file["file1/VDF"] = vcellf
   end
end

main()
```