# Sample postprocessing script for virtual satellite tracking.
#
# Hongyang Zhou, hyzhou@umich.edu 01/25/2021

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
   meta = read_meta(filenames[i])
   id = get_cellid(meta, loc)
   data_series[i] = read_variable_select(meta, var, id)[1][1]
end

# Save into text file
open("satellite.txt", "w") do io
   writedlm(io, data_series)
end

# Specify starting and end point for line tracking
point1 = [12Re, 0, 0]
point2 = [15Re, 0, 0]

# Extract data along the line segment in one snapshot
meta = read_meta(filenames[1])
cellIDs, distances, coords = get_cell_in_line(meta, point1, point2)
