# Sample postprocessing script for virtual satellite tracking.
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, DelimitedFiles, Vlasiator

Re = Vlasiator.Re # Earth radius

# data directory
dir = "./"

filenames = glob("bulk1*.vlsv", dir)

# variable to be extracted
var = "proton/vg_rho"
# virtual satellite location
loc = [12Re, 0, 0]

# Allocate vector for data
data_series = zeros(Float32, length(filenames))

# Extract data from each frame
Threads.@threads for i = 1:length(filenames)
   meta = load(filenames[i])
   id = getcell(meta, loc)
   data_series[i] = readvariable(meta, var, id)[1][1]
end

# Save into text file
open("satellite.txt", "w") do io
   writedlm(io, data_series)
end

# Specify starting and end point for line tracking
point1 = [12Re, 0, 0]
point2 = [15Re, 0, 0]

# Extract data along the line segment in one snapshot
meta = load(filenames[1])
cellIDs, distances, coords = getcellinline(meta, point1, point2)

## Visualization
#=
using PyPlot, DelimitedFiles

data = readdlm("virtual_satellite.txt")

plot(data)
=#