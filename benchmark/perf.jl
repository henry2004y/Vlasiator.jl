using Pkg
Pkg.activate(".")
using Vlasiator
using BenchmarkTools

#file = "bulk.singleprecision.vlsv" # 80 KB rho
#file = "bulk.0000003.vlsv"        # 900 KB rho
file = "bulk1.0001000.vlsv"       # 32 MB rho
#file = "test/data/bulk.2d.vlsv"

@time meta = load(file);

@benchmark meta = load($file)

meta = load(file);

var = "proton/vg_rho";

@benchmark rho = meta[$var]

@benchmark rho_unsorted = readvariable($meta, $var, false)

println("")