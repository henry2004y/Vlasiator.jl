#using Pkg
#Pkg.activate(".")
using Vlasiator
using BenchmarkTools

#filename = "bulk.singleprecision.vlsv"
#filename = "bulk.0000003.vlsv"
filename = "bulk1.0001000.vlsv"
#filename = "test/data/bulk.2d.vlsv"

@time meta = load(filename);

@btime meta = load($filename);

meta = load(filename);

var = "proton/vg_rho";

@btime rho = meta[$var];

@btime rho_unsorted = readvariable($meta, $var, false);