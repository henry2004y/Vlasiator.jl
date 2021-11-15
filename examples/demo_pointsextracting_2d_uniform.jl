# Sample postprocessing script for virtual satellite trackings.
# Simultaneous virtual satellite trackings at uniformly sampled locations on a 2D plane.
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

"Select cells in 2D `meta` with uniform distance `dx`."
function sample(meta, dx)
   dcell = floor(Int, dx ÷ meta.dcoord[1])

   cellid = meta.cellid[meta.cellindex]
   cellid = reshape(cellid, meta.ncells[1], meta.ncells[2])
   cellid_select = cellid[1+dcell:dcell:end-dcell, 1+dcell:dcell:end-dcell]
end

function extract_vars(files, dx)
   nfiles = length(files)

   meta = load(files[1])
   ids = sample(meta, dx)
   locations = [getcellcoordinates(meta, id) for id in ids]

   nsize = size(ids)
   println("Number of virtual satellites: $(length(ids))")

   # variables to be extracted
   t   = zeros(Float32, nfiles)
   rho = zeros(Float32, nfiles, nsize[1], nsize[2])
   vx  = zeros(Float32, nfiles, nsize[1], nsize[2])
   vy  = zeros(Float32, nfiles, nsize[1], nsize[2])
   p   = zeros(Float32, nfiles, nsize[1], nsize[2])
   bz  = zeros(Float32, nfiles, nsize[1], nsize[2])
   ex  = zeros(Float32, nfiles, nsize[1], nsize[2])
   ey  = zeros(Float32, nfiles, nsize[1], nsize[2])

   # Extract data from each frame
   Threads.@threads for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[i,:,:] = readvariable(meta, "proton/vg_rho", ids)
      v = readvariable(meta, "proton/vg_v", ids)
      vx[i,:,:] = v[1,:]
      vy[i,:,:] = v[2,:]
      p[i,:,:] = readvariable(meta, "vg_pressure", ids)
      bz[i,:,:] = readvariable(meta, "vg_b_vol", ids)[3,:]
      e = readvariable(meta, "vg_e_vol", ids)
      ex[i,:,:] = e[1,:]
      ey[i,:,:] = e[2,:]
   end

   # Save into binary file
   jldsave("satellites_uniform_sampled.jld2"; locations, t, rho, vx, vy, p, bz, ex, ey)
end

#####

Re = Vlasiator.Re # Earth radius

files = glob("bulk*.vlsv", "./")

dx = 2Re # uniform sampling distance [m]

println("Number of files: $(length(files))")
println("Running with $(Threads.nthreads()) threads...")

@time extract_vars(files, dx)

println("Virtual satellite extraction done!")