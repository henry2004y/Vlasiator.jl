# Extract variables alone a line across multiple frames.
#
# Usage:
#   julia -t 4 demo_lineextracting_mt.jl
# or
#   JULIA_NUM_THREADS=4 julia demo_lineextracting_mt.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob
using JLD2: jldsave

files = glob("bulk*.vlsv", ".")
nfile = length(files)
# Define end points of the line in Earth radius
x1, x2 = 6.0, 20.0
point1 = [x1, 0, 0] .* Vlasiator.RE
point2 = [x2, 0, 0] .* Vlasiator.RE
# Extract cell info along the line
cellids, distances, coords =
   let meta = load(files[1])
      getcellinline(meta, point1, point2)
   end

n = zeros(Float32, length(cellids), nfile)
v = zeros(Float32, 3, length(cellids), nfile)
p = zeros(Float32, 6, length(cellids), nfile)
b = similar(v)
e = similar(v)

Threads.@threads for i in eachindex(files)
   println("i = $i/$nfile, file = $(files[i])")
   local meta = load(files[i])

   n[:,i] = readvariable(meta, "proton/vg_rho", cellids)
   v[:,:,i] = readvariable(meta, "proton/vg_v", cellids)
   p[1:3,:,i] = readvariable(meta, "proton/vg_ptensor_diagonal", cellids)
   p[4:6,:,i] = readvariable(meta, "proton/vg_ptensor_offdiagonal", cellids)
   b[:,:,i] = readvariable(meta, "vg_b_vol", cellids)
   e[:,:,i] = readvariable(meta, "vg_e_vol", cellids)
end

# Save into binary file
jldsave("vars_sun_earth_line.jld2"; x1, x2, t, n, v, p, b, e)

println("Line extraction finished!")