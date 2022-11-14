# Save Epar, Eperp, B in a selected slice region from 3D VLSV outputs.
#
# Hongyang Zhou, hyzhou@umich.edu 10/31/2022

using Vlasiator, Glob, PyPlot
using JLD2: jldsave
using StaticArrays
using LinearAlgebra: norm

function getvars(meta::MetaVLSV, pArgs::Vlasiator.PlotArgs, normal::Symbol, range1, range2)
   (;origin, idlist, indexlist) = pArgs

   B = readvariable(meta, "vg_b_vol")

   Bout = @views refineslice(meta, idlist, B[:,indexlist], normal)[:, range1, range2]

   ncells = meta.ncells .* 2^meta.maxamr
   if normal == :x
      dir = 1
   elseif normal == :y
      dir = 2
   elseif normal == :z
      dir = 3
   else
      @error "Unknown normal direction $normal"
   end

   sliceratio = (origin - meta.coordmin[dir]) / (meta.coordmax[dir] - meta.coordmin[dir])
   0.0 ≤ sliceratio ≤ 1.0 || error("slice plane index out of bound!")
   # Find the cut plane index for each refinement level
   icut = floor(Int, sliceratio*ncells[dir]) + 1

   E = readvariable(meta, "fg_e", true, true)

   E = if normal == :x
      E[:,icut,range1,range2]
   elseif normal == :y
      E[:,range1,icut,range2]
   elseif normal == :z
      E[:,range1,range2,icut]
   end

   n2D = meta["n"][indexlist]

   n = refineslice(meta, idlist, n2D, normal)[range1, range2]

   v2D = meta["proton/vg_v"][:,indexlist]
   v = @views refineslice(meta, idlist, v2D, normal)[:, range1, range2]

   p2D = meta["P"][indexlist]
   p = refineslice(meta, idlist, p2D, normal)[range1, range2]

   E, Bout, n, v, p
end

function main()
   # Parameters
   directory = "/home/hongyang/runs/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"
   outdir = "EM/"
   outfile = "EM.jld2"
   
   files = glob("bulk1.*vlsv", directory)
   # Only orthogonal slices are supported
   normal = :y # (:x, :y, :z)
   origin = 0.0
   var = "e"
   axisunit = EARTH
   
   meta = load(files[1])
   
   pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
   
   x1, x2 = Vlasiator.get_axis(pArgs)
   
   extent = [-23., 12., -10., 10.]
   
   range1, range2 =
      searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
      searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])

   E = zeros(Float32, 3, length(range1), length(range2), length(files))
   B = zeros(Float32, 3, length(range1), length(range2), length(files))
   n = zeros(Float32, length(range1), length(range2), length(files))
   v = zeros(Float32, 3, length(range1), length(range2), length(files))
   p = zeros(Float32, length(range1), length(range2), length(files))

   for i in eachindex(files)
      meta = load(files[i])
      E[:,:,:,i], B[:,:,:,i], n[:,:,i], v[:,:,:,i], p[:,:,i] =
         getvars(meta, pArgs, normal, range1, range2)
   end

   # Save into binary file
   jldsave(outdir*outfile; E, B, n, v, p)
end

main()