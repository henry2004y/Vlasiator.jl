# Sample script for plotting the phase space density near a given spatial
# location.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

filename = "bulk.0000001.vlsv"

meta = readmeta(filename)

coordinates = [0.0, 0.0, 0.0]

plot_vdf(meta, coordinates; verbose=true)

# Show the spatial distribution of cells with saved VDF
cellswithVDF = getcellwithvdf(meta)
cellID = readvariable(meta, "CellID")
cellID = reshape(cellID, meta.ncells[1], meta.ncells[2])
f_saved_index_ = [findfirst(==(cid), cellID) for cid in cellswithVDF]

x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1]) ./ Vlasiator.Re
y = LinRange(meta.coordmin[2], meta.coordmax[2], meta.ncells[2]) ./ Vlasiator.Re

xcell = zeros(size(cellswithVDF))
ycell = similar(xcell)

for i in eachindex(f_saved_index_)
   xcell[i] = x[f_saved_index_[i][1]]
   ycell[i] = y[f_saved_index_[i][2]]
end

figure()
pcolormesh(meta, "vg_pressure", colorscale=Linear)
ax = plt.gca()
ax.scatter(xcell, ycell, marker="+", color="w")