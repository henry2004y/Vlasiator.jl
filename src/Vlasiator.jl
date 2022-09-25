module Vlasiator

using Requires
using StaticArrays: SVector, @SVector, SMatrix, @SMatrix
using Printf: @sprintf
using LinearAlgebra: ×, dot, ⋅, norm, normalize, normalize!
using Statistics: mean
using EzXML
using Mmap: mmap
using WriteVTK
using LazyGrids: ndgrid
using LaTeXStrings
using Dates
using SnoopPrecompile

include("utility/rotation.jl")
include("utility/log.jl")
include("utility/monitor.jl")
include("utility/differential.jl")
include("utility/fluxfunction.jl")
include("vlsv/vlsvreader.jl")
include("vlsv/vlsvutility.jl")
include("utility/plot.jl")

export
   # vlsvreader
   MetaVLSV, VarInfo,
   load, readvariable, readparameter, readvariablemeta, readvcells, extractsat,
   hasvariable, hasparameter, hasname,
   # vlsvutility
   getcell, getslicecell, getlevel, refineslice, getcellcoordinates, getvcellcoordinates,
   getchildren, getparent, isparent, getsiblings,
   getcellinline, getnearestcellwithvdf, getcellwithvdf,
   getdensity, getvelocity, getpressure, getmaxwellianity, write_vtk, write_vlsv, issame,
   # plot helper
   SI, EARTH, Log, Linear, SymLog,
   # log
   readlog,
   # physical parameter monitor
   check_plasma_characteristics,
   # fluxfunction
   compute_flux_function, find_reconnection_points

# SnoopPrecompile
@precompile_setup begin
   initfile = joinpath(@__DIR__, "../test/init.vlsv")
   @precompile_all_calls begin
      meta = load(initfile)
      readvariable(meta, "CellID")
      readvariable(meta, "proton/vg_v")
      readvariable(meta, "fg_b")
      readvariable(meta, "CellID", UInt(1))
      readvariable(meta, "CellID", UInt[1])
      getcell(meta, [1.0, 1.0, 1.0])
   end
end

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      include("plot/plots.jl")
   end
end

end