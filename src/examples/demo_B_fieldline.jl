# Sample script for plotting fieldlines with handpicked seeds.
#
# Hongyang Zhou, hyzhou@umich.edu

# FieldTracer is not a registered pkg yet
# using Pkg
# Pkg.add(url="https://github.com/henry2004y/FieldTracer.jl.git", rev="master")

using Vlasiator, PyPlot, FieldTracer

filename = "bulk.0000999.vlsv"
meta = readmeta(filename)
plot_pcolormesh(meta, "proton/vg_rho", islinear=true)

Re = Vlasiator.Re
# regular Cartesian mesh
gridx = range(meta.xmin, meta.xmax, length=meta.xcells) 
gridy = range(meta.zmin, meta.zmax, length=meta.zcells)

ns1 = 14 # number of seeding points on the dayside
ns2 = 10 # number of seeding points on the nightside
ns4 = 4  # number of seeding points in the polar region
seeds = Matrix{Float64}(undef, 2, ns1+2*ns2+2*ns4)
for i in 1:ns1
   seeds[1,i] = meta.xmax / ns1 * (i - 1)
   seeds[2,i] = 0.0
end
for i in 1:ns2
   seeds[1,ns1+i] = -10Re
   seeds[2,ns1+i] = meta.zmin + (meta.zmax - meta.zmin) / ns2 * (i - 1)
   seeds[1,ns1+ns2+i] = -30Re
   seeds[2,ns1+ns2+i] = meta.zmin + (meta.zmax - meta.zmin) / ns2 * (i - 1)
end
for i in 1:ns4
   seeds[1,ns1+2*ns2+i] = -20Re + 10Re*(i - 1)
   seeds[2,ns1+2*ns2+i] = 30Re
   seeds[1,ns1+2*ns2+ns4+i] = -20Re + 10Re*(i - 1)
   seeds[2,ns1+2*ns2+ns4+i] = -30Re
end

b = readvariable(meta, "vg_b_vol")
bx = reshape(b[1,:], meta.xcells, meta.zcells)
bz = reshape(b[3,:], meta.xcells, meta.zcells)

for i = axes(seeds,2)
   startx, starty = seeds[:,i]
   x1, y1 = trace2d(bx, bz, startx, starty, gridx, gridy;
      ds=0.5, maxstep=3000, gridType="ndgrid")
   x1 ./= Re
   y1 ./= Re
   if length(x1) < 5; continue; end
   line = plot(x1, y1, color="w")
   add_arrow(line[1])
end