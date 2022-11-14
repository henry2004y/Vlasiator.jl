# Finding X-points in a 2D magnetic reconnection configuration and saving the coordinates
# as well as extracted reconnection rate Ey from multiple outputs.
# This example assumes X-Z meridional plane.
#
# Hongyang Zhou, hyzhou@umich.edu

using JLD2: jldsave
using Vlasiator, Glob, ProgressMeter

function main()
   files = glob("bulk*.vlsv", ".")

   nG = 2 # number of ghost cells

   x, z, dx = load(files[1]) do meta
      LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1]) ./ Vlasiator.RE,
      LinRange(meta.coordmin[3], meta.coordmax[3], meta.ncells[3]) ./ Vlasiator.RE,
      [meta.dcoord[1], meta.dcoord[3]]
   end

   xmin_ = searchsortedfirst(x, 7.0)
   xmax_ = searchsortedlast(x, 9.0)
   zmin_ = searchsortedfirst(z, -4.0)
   zmax_ = searchsortedlast(z, 4.0)

   x_points_x = Vector{Vector{eltype(x)}}(undef, 0)
   x_points_z = similar(x_points_x)

   @showprogress 5 "Finding X-points..." for ifile in eachindex(files)
      meta = load(files[ifile])
      b = meta["vg_b_vol"]
      b = reshape(b, 3, meta.ncells[1], meta.ncells[3])

      flux = compute_flux_function(b, dx, nG)
      indices_x, _ = find_reconnection_points(flux[xmin_:xmax_,zmin_:zmax_], 5e-3)

      push!(x_points_x, x[indices_x[1,:].+xmin_.-1])
      push!(x_points_z, z[indices_x[2,:].+zmin_.-1])
   end

   ## Extract Ey at X-points
   ey = Vector{Vector{Float32}}(undef, 0)

   for it in eachindex(x_points_x)
      meta = load(files[it])
      ids = Vector{Int}(undef, length(x_points_x[it]))
      for ip in eachindex(x_points_x[it])
         loc = [x_points_x[it][ip], 0.0, x_points_z[it][ip]] .* Vlasiator.RE
         ids[ip] = getcell(meta, loc)
      end
      ey_now = readvariable(meta, "vg_e_vol", ids)[2,:]
      push!(ey, ey_now)
   end

   # save
   jldsave("x_point_locations.jld2"; x_points_x, x_points_z, ey)
end

main()