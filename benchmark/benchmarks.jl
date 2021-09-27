using BenchmarkTools, LazyArtifacts

t = @elapsed using Vlasiator
println("Julia version is $VERSION")
println(string("Vlasiator.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Vlasiator.jl...")
println()

directory = artifact"testdata"
files = ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv")

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["variable"])
file = joinpath(directory, files[2])
meta = load(file)
SUITE["read"]["meta"] = @benchmarkable load($file)
SUITE["read"]["DCCRG"] = @benchmarkable readvariable($meta, "proton/vg_rho")
ids = 3000:6300
SUITE["read"]["DCCRG_select"] = @benchmarkable readvariable($meta, "proton/vg_rho", $ids)
SUITE["read"]["FG"] = @benchmarkable readvariable($meta, "fg_b")
file = joinpath(directory, files[1])
meta = load(file)
SUITE["read"]["VDF"] = @benchmarkable readvcells($meta, 2; pop="proton")

SUITE["VTK"] = BenchmarkGroup(["conversion"])
file = joinpath(directory, files[3])
meta = load(file)
SUITE["VTK"]["AMR"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"])