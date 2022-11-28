using Vlasiator
using BenchmarkTools
using Glob

## Reading

#file = "bulk.singleprecision.vlsv" # 80 B rho
#file = "bulk.0000003.vlsv"        # 900 KB rho
file = "bulk1.0001000.vlsv"       # 32 MB rho
#file = "test/data/bulk.2d.vlsv"

@time meta = load(file);
@time meta = load(file)

meta = load(file);

var = "proton/vg_rho";

@time rho = meta[var];
@time rho = meta[var];

@time rho_unsorted = readvariable($meta, $var, false);
@time rho_unsorted = readvariable($meta, $var, false);

## Static location extracting

dir = "/wrk/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"

filenames = glob("bulk1*.vlsv", dir)

var = "proton/vg_rho"
loc = [8Vlasiator.RE, 0, 0]

id = let
   meta = load(filenames[1])
   getcell(meta, loc)
end

@time data_series = extractsat(filenames, var, id);

## Plotting with PyPlot

# 2D scalar contour, uniform mesh
filename = "bulk.0000501.vlsv"
meta = load(filename)
@time pcolormesh(meta, "rho", colorscale=Log)
close()
@time pcolormesh(meta, "rho", colorscale=Log)

# 2D scalar contour, AMR mesh
filename = "bulk1.0001000.vlsv"
meta = load(filename)
@time pcolormesh(meta, "proton/vg_rho", colorscale=Log)
close()
@time pcolormesh(meta, "proton/vg_rho", colorscale=Log)

## Virtual satellite tracking

using Vlasiator: RE
# EGI data on Turso, UH
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