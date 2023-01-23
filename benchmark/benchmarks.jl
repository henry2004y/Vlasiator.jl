using BenchmarkTools, LazyArtifacts

t = @elapsed using Vlasiator
println("Julia version is $VERSION")
println(string("Vlasiator.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Vlasiator.jl...")
println()

#directory = artifact"testdata" # The artifact"" macro does not support a path input!
artifact_toml = LazyArtifacts.find_artifacts_toml(".")
directory = artifact_hash("testdata", artifact_toml) |> artifact_path

files = ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv")

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["variable"])
file = joinpath(directory, files[2])
meta = load(file)
SUITE["read"]["meta"] = @benchmarkable load($file)
SUITE["read"]["DCCRG_scalar"] = @benchmarkable readvariable($meta, "proton/vg_rho")
SUITE["read"]["DCCRG_vector"] = @benchmarkable readvariable($meta, "proton/vg_v")
SUITE["read"]["DCCRG_derived"] = @benchmarkable readvariable($meta, "Vperp")
ids = collect(100:110)
SUITE["read"]["DCCRG_select_small"] =
   @benchmarkable readvariable($meta, "proton/vg_rho", $ids)
ids = collect(3000:6300)
SUITE["read"]["DCCRG_select_large"] =
   @benchmarkable readvariable($meta, "proton/vg_rho", $ids)
SUITE["read"]["FG"] = @benchmarkable readvariable($meta, "fg_b")
file = joinpath(directory, files[1])
meta = load(file)
SUITE["read"]["VDF"] = @benchmarkable readvcells($meta, 5; species="proton")

SUITE["VTK"] = BenchmarkGroup(["conversion"])
file = joinpath(directory, files[3])
meta = load(file)
SUITE["VTK"]["AMR_full"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"])
SUITE["VTK"]["AMR_maxonly"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"],
   skipghosttype=true, maxamronly=true)