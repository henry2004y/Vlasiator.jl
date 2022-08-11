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

precompile(load, (String,))
precompile(readvariable, (MetaVLSV, String, Vector{UInt}))
precompile(readvariable, (MetaVLSV, String, Int))
precompile(readvariable, (MetaVLSV, String))
precompile(readvariable, (MetaVLSV, String, UnitRange{Int}))
precompile(readvariable, (MetaVLSV, String, Bool))
precompile(readvcells, (MetaVLSV, Int))
precompile(readvcells, (MetaVLSV, UInt))
precompile(readlog, (String,))
precompile(readmesh, (IOStream, EzXML.Node, String, String))
precompile(readparameter, (MetaVLSV, String))
precompile(readparameter, (IOStream, EzXML.Node, String))
precompile(readvariablemeta, (MetaVLSV, String))
precompile(readvector, (IOStream, EzXML.Node, String, String))
precompile(getcell, (MetaVLSV, Vector{Float64}))
precompile(getObjInfo, (EzXML.Node, String, String, String))
precompile(get_axis, (Vlasiator.PlotArgs,))
precompile(getcellwithvdf, (MetaVLSV,))
precompile(getdata2d, (MetaVLSV, String))
precompile(getfooter, (IOStream,))
precompile(hasname, (EzXML.Node, String, String))
precompile(hasparameter, (MetaVLSV, String))
precompile(hasvariable, (MetaVLSV, String))
precompile(prep2d, (MetaVLSV, String))
precompile(set_args, (MetaVLSV, String, Vlasiator.AxisUnit))
precompile(searchsorted, (Vector{UInt}, Int))

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      include("plot/plots.jl")
   end
end

end