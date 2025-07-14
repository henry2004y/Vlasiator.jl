using BenchmarkTools, Downloads, Tar, CodecZlib

println("Julia version is $VERSION")
println()
println("Benchmarking Vlasiator.jl...")
println()

testdata_url = "https://github.com/henry2004y/vlsv_data/raw/master/testdata.tar.gz"
datadir = "data"
files = ("bulk.1d.vlsv", "bulk.2d.vlsv", "bulk.amr.vlsv")

# Check if all files already exist
if joinpath(datadir, files[1]) -> isfile
   println("✅ All data files already exist.")
else
   println("⬇️ Downloading and extracting data...")

   # Download and extract the data
   testdata = Downloads.download(testdata_url)
   open(Gzip.GzipDecompressorStream, testdata) do io
      Tar.extract(io, datadir)
   end

   println("📦 Extraction complete.")
end

const SUITE = BenchmarkGroup()

SUITE["read"] = BenchmarkGroup(["variable"])
file = joinpath(datadir, files[2])
meta = load(file)
SUITE["read"]["meta"] = @benchmarkable load($file)
SUITE["read"]["DCCRG_scalar"] = @benchmarkable readvariable($meta, "proton/vg_rho")
SUITE["read"]["DCCRG_vector"] = @benchmarkable readvariable($meta, "proton/vg_v")
SUITE["read"]["DCCRG_derived"] = @benchmarkable readvariable($meta, "Vperp")
ids = collect(100:110)
SUITE["read"]["DCCRG_select_small"] = @benchmarkable readvariable(
   $meta, "proton/vg_rho", $ids)
ids = collect(3000:6300)
SUITE["read"]["DCCRG_select_large"] = @benchmarkable readvariable(
   $meta, "proton/vg_rho", $ids)
SUITE["read"]["FG"] = @benchmarkable readvariable($meta, "fg_b")
file = joinpath(datadir, files[1])
meta = load(file)
SUITE["read"]["VDF"] = @benchmarkable readvcells($meta, 5; species = "proton")

SUITE["VTK"] = BenchmarkGroup(["conversion"])
file = joinpath(datadir, files[3])
meta = load(file)
SUITE["VTK"]["AMR_full"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"])
SUITE["VTK"]["AMR_maxonly"] = @benchmarkable Vlasiator.fillmesh($meta,
   $["proton/vg_rho", "proton/vg_v", "fg_b", "vg_boundarytype"],
   skipghosttype = true, maxamronly = true)
