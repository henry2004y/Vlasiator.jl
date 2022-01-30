# Sample postprocessing script for virtual satellite trackings.
# Simultaneous virtual satellite trackings at uniformly sampled locations on a 2D plane.
# Outputs are stored in binary format for sharing within Julia.
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
   vz  = zeros(Float32, nfiles, nsize[1], nsize[2])
   p   = zeros(Float32, nfiles, nsize[1], nsize[2])
   bx  = zeros(Float32, nfiles, nsize[1], nsize[2])
   by  = zeros(Float32, nfiles, nsize[1], nsize[2])
   bz  = zeros(Float32, nfiles, nsize[1], nsize[2])
   ex  = zeros(Float32, nfiles, nsize[1], nsize[2])
   ey  = zeros(Float32, nfiles, nsize[1], nsize[2])
   ez  = zeros(Float32, nfiles, nsize[1], nsize[2])

   # Extract data from each frame
   for i = eachindex(files)
      meta = load(files[i])
      t[i] = meta.time
      rho[i,:,:] = readvariable(meta, "proton/vg_rho", ids)
      v = readvariable(meta, "proton/vg_v", ids)
      vx[i,:,:] = v[1,:]
      vy[i,:,:] = v[2,:]
      vz[i,:,:] = v[3,:]
      p[i,:,:] = readvariable(meta, "vg_pressure", ids)
      b = readvariable(meta, "vg_b_vol", ids)
      bx[i,:,:] = b[1,:]
      by[i,:,:] = b[2,:]
      bz[i,:,:] = b[3,:]
      e = readvariable(meta, "vg_e_vol", ids)
      ex[i,:,:] = e[1,:]
      ey[i,:,:] = e[2,:]
      ez[i,:,:] = e[3,:]
   end

   # Save into binary file
   jldsave("satellites_uniform_sampled.jld2";
      locations, t, rho, vx, vy, vz, p, bx, by, bz, ex, ey, ez)
end

#####

files = glob("bulk*.vlsv", "./")

dx = 5Re # uniform sampling distance [m]

println("Number of files: $(length(files))")
println("Running with $(Threads.nthreads()) threads...")

@time extract_vars(files, dx)

println("Virtual satellite extraction done!")