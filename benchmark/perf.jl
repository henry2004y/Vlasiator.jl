# Script for generating benchmark results in
# https://henry2004y.github.io/Vlasiator.jl/dev/benchmark/

using Vlasiator, PyPlot
using BenchmarkTools

# Turn off GUI for Matplotlib
pygui(false)

files = ["1d_single.vlsv", "bulk.2d.vlsv", "2d_double.vlsv", "2d_AFC.vlsv", "3d_EGI.vlsv"]

# Download test files if not found in the current path
for i in eachindex(files)
   if isfile(files[i])
      @info "Benchmark file $(files[i]) found..."
      continue
   elseif i == 4
      url_base = "https://a3s.fi/swift/v1/AUTH_81f1cd490d494224880ea77e4f98490d/vlasiator-2d-afc/"
      filename = "production_halfres/bulk.0000000.vlsv"
      file_origin = basename(filename)
      url = joinpath(url_base, filename)
      @info "Downloading 1.8G test file..."
      run(`curl -o $(files[i]) $url`)
   elseif i in (1,2,3)
      @info "Downloading test files..."
      run(`curl -o testdata.tar.gz https://raw.githubusercontent.com/henry2004y/vlsv_data/master/testdata.tar.gz`)
      run(`curl -o 1d_single.vlsv https://raw.githubusercontent.com/henry2004y/vlsv_data/master/1d_single.vlsv`)
      run(`curl -o 2d_double.vlsv https://raw.githubusercontent.com/henry2004y/vlsv_data/master/2d_double.vlsv`)
      run(`tar -xzvf testdata.tar.gz`)
   elseif i == 5
      @warn "$(files[i]) is not open-access!"
   end
end

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

# Add children groups and tags to our benchmark suite.
suite["load"] = BenchmarkGroup(files)
suite["read"] = BenchmarkGroup(files)
suite["plot"] = BenchmarkGroup(["2d", "3d"])

# Add benchmarks
for (i, file) in enumerate(files)
   !isfile(file) && continue
   suite["load"][file] = @benchmarkable load($file)
   if i == 4
      var = "proton/rho"
      suite["read"][file*"_sorted"] = @benchmarkable meta[$var] setup=(meta=load($file))
      suite["read"][file*"_unsorted"] =
         @benchmarkable readvariable(meta, $var, false) setup=(meta=load($file))
      suite["plot"]["contour_uniform"] =
         @benchmarkable pcolormesh(meta, $var, colorscale=Log) setup=(meta=load($file))
   else
      var = "proton/vg_rho"
      suite["read"][file*"_sorted"] = @benchmarkable meta[$var] setup=(meta=load($file))
      suite["read"][file*"_unsorted"] =
         @benchmarkable readvariable(meta, $var, false) setup=(meta=load($file))
   end
   if i == 5
      var = "proton/vg_rho"
      suite["plot"]["contour_nonuniform"] =
         @benchmarkable pcolormesh(meta, $var, colorscale=Log) setup=(meta=load($file))
   end
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
   loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals)
else
   tune!(suite)
   BenchmarkTools.save(paramspath, params(suite))
end

results = run(suite, verbose=true, samples=100, seconds=20)

show(results)

if isfile(files[5])
   println("")
   println("----------")
   meta = load(files[5])
   # mmap
   @time b = readvariable(meta, "fg_b", true, true);
   println("----------")
   # regular
   @time b = readvariable(meta, "fg_b", true, false);
end
