# Sample script for stream tracing from a given starting point.
#
# Hongyang Zhou, hyzhou@umich.edu 01/27/2021

# FieldTracer is not a registered pkg yet
# using Pkg
# Pkg.add(url="https://github.com/henry2004y/FieldTracer.jl.git", rev="master")
using PyPlot, FieldTracer, Vlasiator

Re = Vlasiator.Re

filename = "bulk.0000501.vlsv"
nameρ = "rho"
nameV = "rho_v"

meta = read_meta(filename)

plot_pcolormesh(meta, nameρ)

v = read_variable(meta, nameV)
vx = reshape(v[1,:], meta.xcells, meta.ycells)
vy = reshape(v[2,:], meta.xcells, meta.ycells)
# tracing starting point
xstart, ystart = 12Re, 0Re
# regular Cartesian mesh
x = range(meta.xmin, meta.xmax, length=meta.xcells) 
y = range(meta.ymin, meta.ymax, length=meta.ycells)

# RK4 scheme by default
x1, y1 = trace2d(vx, vy, xstart, ystart, x, y;
   ds=0.5, maxstep=3000, gridType="ndgrid")
x1 ./= Re
y1 ./= Re

plot(x1, y1)

streamline(meta, nameV, comp="xy", color="w", density=1.0)