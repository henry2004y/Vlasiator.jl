# Sample script for plotting the phase space density near a given spatial location.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

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
   xcell[i] = locations[i][1] / Vlasiator.Re
   ycell[i] = locations[i][2] / Vlasiator.Re
end

figure()
pcolormesh(meta, "vg_pressure", colorscale=Linear)
ax = plt.gca()
ax.scatter(xcell, ycell, marker="+", color="w")