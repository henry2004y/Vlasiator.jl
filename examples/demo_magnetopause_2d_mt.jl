# Extracting variable values along the magnetopause Bz==0 contour from multiple frames.
#
# Usage:
#   julia -t 4 demo_magnetopause_2d_mt.jl
# or
#   JULIA_NUM_THREADS=4 julia demo_magnetopause_2d_mt.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob
using Vlasiator: RE # Earth radius, [m]
using JLD2: jldsave

"""
    extract_magnetopause_var(files; zmin=-4RE, zmax=4RE, verbose=true)

Extract variables on the magnetopause defined by ``B_z = 0``.
"""
function extract_magnetopause_var(files; zmin=-4RE, zmax=4RE, verbose=true)
   nfiles = length(files)

   verbose && println("Number of files: $nfiles")

   meta = load(files[1])

   x = LinRange{Float32}(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])
   z = LinRange{Float32}(meta.coordmin[3], meta.coordmax[3], meta.ncells[3])

   z_range = let zmin_ = searchsortedfirst(z, zmin), zmax_ = searchsortedlast(z, zmax)
      zmin_:zmax_
   end

   x_crossing = zeros(Float32, length(z_range))
   z_crossing = z[z_range]

   cellids = Vector{Int}(undef, length(z_range))

   x_crossings = zeros(length(z_range), nfiles)
   v = zeros(3, length(z_range), nfiles)
   ey = zeros(length(z_range), nfiles)

   Threads.@threads for ifile in eachindex(files)
      verbose && println("$(files[ifile]) on thread $(Threads.threadid())")
      meta = load(files[ifile])

      # Obtain thermal temperature
      b = meta["vg_b_vol"]
      b = reshape(b, 3, meta.ncells[1], meta.ncells[3])

      # Extract the last point from right to left which fulfills Bz < 0
      for (i,k) in enumerate(z_range) # scan in z direction
         ind_ = findlast(>(0), @view b[3,:,k]) + 1 # count from upstream
         isnothing(ind_) && (ind_ = 1) # if not found then set to 1
         x_crossing[i] = x[ind_]
      end

      cellids = [ getcell(meta, [x_crossing[i], 0, z_crossing[i]])
         for i in eachindex(x_crossing, z_crossing) ]

      x_crossings[:,ifile] = x_crossing
      v[:,:,ifile] = readvariable(meta, "proton/vg_v", cellids)
      ey[:,ifile] = readvariable(meta, "vg_e_vol", cellids)[2,:]
   end

   z_crossing, x_crossings, v, ey
end

######

files = glob("bulk*.vlsv", ".")

@time z,x,v,ey = extract_magnetopause_var(files)

jldsave("magnetopause.jld2"; z,x,v,ey)