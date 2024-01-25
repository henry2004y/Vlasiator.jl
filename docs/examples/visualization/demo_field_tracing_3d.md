# ---
# title: Tracing field lines in 3D
# id: demo_3d_trace
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to trace streamlines in a 3D field
# ---

This demo shows how to trace along a field line and extract variables from 3D AMR grid.
!!! warning
    This is a proof of concept, which is memory-bound in 3D and may be quite inefficient!

```julia
using Vlasiator, FieldTracer, PyPlot
using StaticArrays
using Vlasiator: RE

function get_cell_along_B(meta::MetaVLSV, seed::Vector{T};
   ds::Float64=0.2, maxstep::Int=20000, direction::String="both") where T
   B = Vlasiator.fillmesh(meta, ["vg_b_vol"]; maxamronly=true)[1][1][1]

   nx, ny, nz = size(B)[2:4]

   xrange = LinRange(meta.coordmin[1], meta.coordmax[1], nx)
   yrange = LinRange(meta.coordmin[2], meta.coordmax[2], ny)
   zrange = LinRange(meta.coordmin[3], meta.coordmax[3], nz)

   bx = @view B[1,:,:,:]
   by = @view B[2,:,:,:]
   bz = @view B[3,:,:,:]
   xs, ys, zs = seed

   x, y, z = trace(bx, by, bz, xs, ys, zs, xrange, yrange, zrange; ds, maxstep, direction)

   cellids = [getcell(meta, SVector(x[1], y[1], z[1]))]

   for i in eachindex(x)[2:end]
      cellidnew = getcell(meta, SVector(x[i], y[i], z[i]))
      if cellidnew != cellids[end]
         push!(cellids, cellidnew)
      end
   end

   cellids
end

function main()
   file = "bulk1.0001000.vlsv"
   meta = load(file)

   seed = [-5.0, 0., 3.5] .* RE
   # For EGI, ds=1.0 is between 1-2 steps per finest cell.
   cellids = get_cell_along_B(meta, seed; ds=1.0, maxstep=100, direction="backward")

   #n = readvariable(meta, "proton/vg_rho", cellids)
   #v = readvariable(meta, "proton/vg_v", cellids)
   #E = readvariable(meta, "fg_e", cellids) # Not available for now!

   x = zeros(Float64, length(cellids))
   y = zeros(Float64, length(cellids))
   z = zeros(Float64, length(cellids))

   for i in eachindex(cellids)
      coords = getcellcoordinates(meta, cellids[i])
      x[i], y[i], z[i] = coords
   end

   fig = plt.figure()
   ax = fig.add_subplot(projection="3d")

   ax.set_xlabel("x [RE]")
   ax.set_ylabel("y [RE]")
   ax.set_zlabel("z [RE]")

   ax.set_box_aspect([1,1,1])
   ax.set_xlim3d([-20.0, -5.0])
   ax.set_ylim3d([0.0, 0.5])
   ax.set_zlim3d([-2.0, 3.0])

   ax.plot3D(x./RE, y./RE, z./RE)
end

main()
```