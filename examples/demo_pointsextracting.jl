# Sample postprocessing script for virtual satellite trackings.
# Simultaneous virtual satellite trackings at multiple locations.
# Outputs are stored in binary format for sharing within Julia.
#
# Usage:
#   julia -t nthreads demo_pointsextracting.jl
# or
#   JULIA_NUM_THREADS=nthreads julia demo_pointsextracting.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Vlasiator
using Vlasiator: RE # Earth radius, [m]
using JLD2: jldsave

function extract_vars(files, locations)
   nfiles = length(files)
   nlocs = length(locations)
   # variables to be extracted
   t   = zeros(Float32, nfiles)
   rho = zeros(Float32, nlocs, nfiles)
   v   = zeros(Float32, 3, nlocs, nfiles)
   p   = zeros(Float32, nlocs, nfiles)
   b   = zeros(Float32, 3, nlocs, nfiles)
   e   = zeros(Float32, 3, nlocs, nfiles)

   ids = Vector{Int}(undef, nlocs)
   load(files[1]) do meta
      for iloc in eachindex(locations)
         ids[iloc] = getcell(meta, locations[iloc])
      end
   end

   # Extract data from each frame
   Threads.@threads for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[:,i] = readvariable(meta, "proton/vg_rho", ids)
      v[:,:,i] = readvariable(meta, "proton/vg_v",   ids)
      p[  :,i] = readvariable(meta, "vg_pressure",   ids)
      b[:,:,i] = readvariable(meta, "vg_b_vol",      ids)
      e[:,:,i] = readvariable(meta, "vg_e_vol",      ids)
   end

   # Save into binary file
   jldsave("satellites.jld2"; locations, t, rho, v, p, b, e)
end

#####

files = glob("bulk*.vlsv", "./")

# virtual satellite locations
locations = [[7RE, 0, 0], [9RE, 0, 0], [11RE, 0, 0], [12RE, 0, 0], [13RE, 0, 0],
   [14RE, 0, 0], [15RE, 0, 0], [16RE, 0, 0], [17RE, 0, 0], [29.3RE, 0, 0]]

println("Number of files: $(length(files))")
println("Number of virtual satellites: $(length(locations))")
println("Running with $(Threads.nthreads()) threads...")

@time extract_vars(files, locations)

println("Virtual satellite extraction done!")