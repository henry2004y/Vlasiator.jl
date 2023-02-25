# ---
# title: Extract bow shock location
# id: demo_2d_bowshock_mt
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to extract the bow shock location and save into file
# ---

This demo shows how to extract the bow shock location from 2D equatorial run outputs and save into file.
To run in multi-threading mode,
```shell
julia -t nthreads demo_bowshock_2d_mt.jl
```
or alternatively
```shell
JULIA_NUM_THREADS=$nthreads julia demo_bowshock_2d_mt.jl
```

```julia
using Vlasiator, Glob
using Vlasiator: RE # Earth radius, [m]
using JLD2: jldsave

# Upstream solar wind temperature
const Tsw = 0.5e6 #[K]

function extract_bowshock_position(files; verbose=true)
   nfiles = length(files)

   verbose && println("Number of files: $nfiles")

   meta = load(files[1])

   x = LinRange{Float32}(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
   y = LinRange{Float32}(meta.coordmin[2], meta.coordmax[2], meta.ncells[2])

   close(meta.fid)

   # Only extract bow shock location near the front between y = ±20RE
   ymin_ = findlast(<(-20RE), y)
   ymax_ = findfirst(>(20RE), y)

   x_crossing = zeros(Float32, ymax_-ymin_+1, nfiles)
   y_crossing = y[ymin_:ymax_]

   Threads.@threads for ifile = 1:nfiles
      verbose && println("$(files[ifile]) on thread $(Threads.threadid())")
      f = load(files[ifile])
      # Obtain thermal temperature
      T = f["T"]
      close(f.fid)
      T = reshape(T, f.ncells[1], f.ncells[2])
      # Extract bow shock location from the 1st point which fulfills
      # the threshold: T > 4 * Tsw
      for (i,j) in enumerate(ymin_:ymax_) # scan in y direction
         ind_ = findlast(>(4*Tsw), @view T[:,j]) # count from upstream
         x_crossing[i,ifile] = x[ind_]
      end
   end

   return x_crossing, y_crossing
end

#####
files = glob("bulk*.vlsv", ".")

@time x_crossing, y_crossing = extract_bowshock_position(files)

jldsave("example.jld2"; x_crossing, y_crossing)
```