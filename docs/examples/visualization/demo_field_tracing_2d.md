# ---
# title: Tracing field lines in 2D
# id: demo_2d_trace
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to trace streamlines in a 2D field
# ---

This demo shows how to trace streamlines in a 2D uniform field.
```julia
using VlasiatorPyPlot, FieldTracer
using Vlasiator: RE # Earth radius, [m]

function main()
   file = "bulk.0000501.vlsv"
   nameρ = "rho"
   nameV = "rho_v"

   meta = load(file)

   pcolormesh(meta, nameρ)

   v = readvariable(meta, nameV)
   vx = reshape(v[1,:], meta.ncells[1], meta.ncells[2])
   vy = reshape(v[2,:], meta.ncells[1], meta.ncells[2])
   # tracing starting point
   xstart, ystart = 12RE, 0RE
   # regular Cartesian mesh
   x = range(meta.coordmin[1], meta.coordmax[1], length=meta.ncells[1])
   y = range(meta.coordmin[2], meta.coordmax[2], length=meta.ncells[2])

   # RK4 scheme by default
   x1, y1 = trace2d(vx, vy, xstart, ystart, x, y;
      ds=0.5, maxstep=3000, gridType="ndgrid")
   x1 ./= RE
   y1 ./= RE

   plot(x1, y1)

   streamplot(meta, nameV, comp="xy", color="w", density=1.0)
end

main()
```