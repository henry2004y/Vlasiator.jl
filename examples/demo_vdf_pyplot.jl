# Sample script for plotting the phase space density near a given spatial location.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot
using Vlasiator: RE # Earth radius [m]

function main()
   file = "bulk.0001347.vlsv"
   meta = load(file)

   coordinates = [0.0, 0.0, 0.0]

   vdfslice(meta, coordinates; verbose=true)

   # Show the spatial distribution of cells with saved VDF
   cellswithVDF = getcellwithvdf(meta)
   locations = [getcellcoordinates(meta, cid) for cid in cellswithVDF]

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