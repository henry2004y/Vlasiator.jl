# Sample postprocessing script for virtual satellite trackings.
# Simultaneous virtual satellite trackings at multiple locations.
# Outputs are stored in binary format for sharing within Julia.
#
# Usage:
#   julia -t nthreads demo_virtual_satellites.jl
# or
#   JULIA_NUM_THREADS=nthreads julia demo_virtual_satellites.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Glob, Vlasiator
using JLD2: jldsave

function extract_vars(files, locations)
   nfiles = length(files)
   nlocs = length(locations)
   # variables to be extracted
   t   = zeros(Float32, nfiles)
   rho = zeros(Float32, nfiles, nlocs)
   vx  = zeros(Float32, nfiles, nlocs)
   vy  = zeros(Float32, nfiles, nlocs)
   p   = zeros(Float32, nfiles, nlocs)
   bz  = zeros(Float32, nfiles, nlocs)
   ex  = zeros(Float32, nfiles, nlocs)
   ey  = zeros(Float32, nfiles, nlocs)

   ids = Vector{UInt64}(undef, nlocs)
   let meta = load(files[1])
      for iloc in eachindex(locations)
         ids[iloc] = getcell(meta, locations[iloc])
      end
   end

   # Extract data from each frame
   Threads.@threads for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[i,:] = readvariable(meta, "proton/vg_rho", ids)
      v = readvariable(meta, "proton/vg_v", ids)
      vx[i,:] = v[1,:]
      vy[i,:] = v[2,:]
      p[i,:] = readvariable(meta, "vg_pressure", ids)
      bz[i,:] = readvariable(meta, "vg_b_vol", ids)[3,:]
      e = readvariable(meta, "vg_e_vol", ids)
      ex[i,:] = e[1,:]
      ey[i,:] = e[2,:]
   end

   # Save into binary file
   jldsave("satellites.jld2"; t, rho, vx, vy, p, bz, ex, ey)
end

#####

Re = Vlasiator.Re # Earth radius

files = glob("bulk*.vlsv", "./")

# virtual satellite locations
loc = [[7Re, 0, 0], [9Re, 0, 0], [11Re, 0, 0], [12Re, 0, 0], [13Re, 0, 0], [14Re, 0, 0],
   [15Re, 0, 0], [16Re, 0, 0], [17Re, 0, 0], [29.3Re, 0, 0]]

println("Number of files: $(length(files))")
println("Number of virtual satellites: $(length(locations))")
println("Running with $(Threads.nthreads()) threads...")

@time extract_vars(files, loc)

println("Virtual satellite extraction done!")