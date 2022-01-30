# Sample script for stream tracing from a given starting point.
#
# Hongyang Zhou, hyzhou@umich.edu

using PyPlot, FieldTracer, Vlasiator
using Vlasiator: RE # Earth radius, [m]

file = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

meta = load(file)

pcolormesh(meta, nameρ)

v = readvariable(meta, nameV)
vx = reshape(v[1,:], meta.ncells[1], meta.ncells[2])
vy = reshape(v[2,:], meta.ncells[1], meta.ncells[2])
# tracing starting point
xstart, ystart = 12Re, 0Re
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