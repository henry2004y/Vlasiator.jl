# Extract the bow shock location from 2D equatorial run.
# Courtesy of Vertti Tarvus, reimplemented from Python. 
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

# Upstream solar wind temperature
Tsw = 0.5e6 #[K]

filename = "bulk.0001347.vlsv"
meta = load(filename)

x = LinRange{Float32}(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
y = LinRange{Float32}(meta.coordmin[2], meta.coordmax[2], meta.ncells[2])

# Obtain thermal temperature
T = meta["T"]
T = reshape(T, meta.ncells[1], meta.ncells[2])

x_crossing = zeros(Float32, meta.ncells[2])
y_crossing = y

# Extract bow shock location from the 1st point which fulfills the threshold: T > 4 * Tsw
# Note: this won't work if the bow shock does not cover the whole domain y extension! 
for j in 1:meta.ncells[2] # scan in y direction
   ind_ = findlast(>(4*Tsw), @view T[:,j]) # count from upstream
   x_crossing[j] = x[ind_]
end

# Julia is column-major, while Python is row-major
imshow(T', extent=(x[1], x[end], y[1], y[end]), origin="lower", cmap=plt.get_cmap("ocean"))
plot(x_crossing, y_crossing, "r")
axis("scaled")
savefig("bs_temp_test.png")