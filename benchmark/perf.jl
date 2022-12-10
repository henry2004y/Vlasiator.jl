# Script for generating benchmark results in
# https://henry2004y.github.io/Vlasiator.jl/dev/benchmark/
#
# Note: this script is not intended for execution as a whole: each benchmark needs to be
# executed independently after importing all the required packages.

using Vlasiator
using Vlasiator: RE
using BenchmarkTools
using Glob

## Reading DCCRG variables

file1 = "bulk.singleprecision.vlsv" # 80 B
file2 = "bulk.0000003.vlsv"         # 900 KB
file3 = "bulk1.0001000.vlsv"        # 32 MB
# Select file
file = file1
# Load metadata
@time meta = load(file);
@time meta = load(file);
@time meta = load(file);
# Select variable name
var = "proton/vg_rho";
# Read variable, sorted
@time rho = meta[var];
@time rho = meta[var];
# Read variable, unsorted
@time rho_unsorted = readvariable($meta, $var, false);
@time rho_unsorted = readvariable($meta, $var, false);

## Plotting with PyPlot
# Log color scale is used for comparison with Analysator.

# 2D density contour on a uniform mesh
filename = "bulk.0000501.vlsv"
meta = load(filename)
@time pcolormesh(meta, "rho", colorscale=Log)
close()
@time pcolormesh(meta, "rho", colorscale=Log)

# 2D density slices from 3D AMR mesh
filename = "bulk1.0001000.vlsv"
meta = load(filename)
@time pcolormesh(meta, "proton/vg_rho", colorscale=Log)
close()
@time pcolormesh(meta, "proton/vg_rho", colorscale=Log)

## Virtual satellite tracking at a static location

# 3D EGI data on Turso, University of Helsinki
dir = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"

filenames = glob("bulk1*.vlsv", dir)

var = "proton/vg_rho"
loc = [12, 0, 0] .* RE

println("Number of files: ", length(filenames))
# Static cell ID
id = let
   meta = load(filenames[1])
   getcell(meta, loc)
end

@btime data_series = extractsat($filenames, $var, $id)