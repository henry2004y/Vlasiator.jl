using Pkg
Pkg.activate(".")
using Vlasiator
using BenchmarkTools

#filename = "bulk.singleprecision.vlsv" # 80 KB rho
#filename = "bulk.0000003.vlsv"        # 900 KB rho
filename = "bulk1.0001000.vlsv"       # 32 MB rho
#filename = "test/data/bulk.2d.vlsv"

@time meta = load(filename);

@benchmark meta = load($filename)

meta = load(filename);

var = "proton/vg_rho";

@benchmark rho = meta[$var]

@benchmark rho_unsorted = readvariable($meta, $var, false)

println("")