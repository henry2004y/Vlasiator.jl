# ---
# title: VDF
# id: demo_vdf
# cover: ../../src/figures/phase_space_distribution.png
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to plot velocity distribution function
# ---

This demo shows how to plot the phase space density (i.e. velocity distribution function) near a given spatial location.
```julia
using VlasiatorPyPlot
using Vlasiator: RE # Earth radius [m]

function main()
   file = "bulk.0001347.vlsv"
   meta = load(file)
   species = "proton"

   coordinates = [0.0, 0.0, 0.0]

   vdfslice(meta, coordinates; verbose=true)

   # Show the spatial distribution of cells with saved VDF
   init_cellswithVDF!(meta, species)
   locations = [getcellcoordinates(meta, cid) for cid in meta.meshes[species].cellwithVDF]

   xcell = zeros(size(locations))
   ycell = similar(xcell)

   for i in eachindex(locations)
      xcell[i] = locations[i][1] / RE
      ycell[i] = locations[i][2] / RE
   end

   figure()
   pcolormesh(meta, "vg_pressure", colorscale=Linear)
   ax = plt.gca()
   ax.scatter(xcell, ycell, marker="+", color="w")
end

main()
```