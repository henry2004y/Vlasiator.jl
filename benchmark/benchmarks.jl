using BenchmarkTools

t = @elapsed using Vlasiator
println("Julia version is $VERSION")
println(string("Vlasiator.jl loading time: \e[33;1;1m$t\e[m seconds"))
println()
println("Benchmarking Vlasiator.jl...")
println()

directory = "test/data"
filenames = ["bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv"]

using ZipFile
r = ZipFile.Reader(joinpath(directory, "bulk_vlsv.zip"))
for file in r.files
   open(joinpath(directory, file.name), "w") do io
      write(io, read(file, String))
   end
end

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["variable"])
filename = joinpath(directory, "bulk.2d.vlsv")
meta = readmeta(filename)
SUITE["read"]["meta"] = @benchmarkable readmeta($filename)
SUITE["read"]["DCCRG"] = @benchmarkable readvariable($meta, "proton/vg_rho")
SUITE["read"]["FG"] = @benchmarkable readvariable($meta, "fg_b")
filename = joinpath(directory, "bulk.1d.vlsv")
meta = readmeta(filename)
SUITE["read"]["VDF"] = @benchmarkable readvcells($meta, 2; pop="proton")

SUITE["VTK"] = BenchmarkGroup(["conversion"])
filename = joinpath(directory, "bulk.amr.vlsv")
meta = readmeta(filename)
SUITE["VTK"]["AMR"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"])