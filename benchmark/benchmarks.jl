using BenchmarkTools

t = @elapsed using Vlasiator
println("Julia version is $VERSION")
println(string("Vlasiator.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Vlasiator.jl...")
println()

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["variable"])
filename = "test/data/bulk.2d.vlsv"
meta = readmeta(filename)
SUITE["read"]["meta"] = @benchmarkable readmeta($filename)
SUITE["read"]["DCCRG"] = @benchmarkable readvariable($meta, "proton/vg_rho")
SUITE["read"]["FG"] = @benchmarkable readvariable($meta, "fg_b")
filename = "test/data/bulk.1d.vlsv"
meta = readmeta(filename)
SUITE["read"]["VDF"] = @benchmarkable readvcells($meta, 2; pop="proton")

SUITE["VTK"] = BenchmarkGroup(["conversion"])
filename = "test/data/bulk.amr.vlsv"
SUITE["VTK"]["image"] = @benchmarkable write_vtk($filename; vti=true)
SUITE["VTK"]["clean"] = @benchmarkable rm("test/data/bulk.amr.vti", force=true)