module Vlasiator

using StaticArrays: SVector, @SVector, SMatrix, @SMatrix
using Printf: @sprintf
using LinearAlgebra: ×, dot, ⋅, norm, normalize, normalize!
using Statistics: mean
using XML: Node, LazyNode, Element, Document, children, tag, attributes, value
import XML
using Mmap: mmap
using LazilyInitializedFields
using WriteVTK
using LazyGrids: ndgrid
using LaTeXStrings
using Dates
using Parsers
using PrecompileTools: @setup_workload, @compile_workload

include("utility/rotation.jl")
include("utility/log.jl")
include("utility/monitor.jl")
include("utility/differential.jl")
include("utility/fluxfunction.jl")
include("vlsv/vlsvreader.jl")
include("vlsv/vlsvvariables.jl")
include("vlsv/vlsvutility.jl")
include("utility/plot.jl")
include("utility/plotrecipe.jl")
include("utility/viz.jl")

export
   # vlsvreader
   MetaVLSV, VarInfo,
   load, readvariable, readparameter, readvariablemeta, readvcells, extractsat,
   hasvariable, hasparameter, hasname,
   # vlsvutility
   getcell, getslicecell, getlevel, refineslice, getcellcoordinates, getvcellcoordinates,
   getchildren, getparent, isparent, getsiblings, getcellinline, getnearestcellwithvdf,
   getmaxwellianity, getKLdivergence, getdensity, getvelocity, getpressure,
   getheatfluxvector, write_vtk, write_vlsv, issame,
   # plot helper
   SI, EARTH, Log, Linear, SymLog,
   # log
   readlog,
   # physical parameter monitor
   check_plasma_characteristics,
   # fluxfunction
   compute_flux_function, find_reconnection_points,
   # visualization
   viz, viz!, vlheatmap, vlslice, vlslices, vdfvolume, vdfslice, vdfslices

# PrecompileTools
@setup_workload begin
   initfile = joinpath(@__DIR__, "../test/init.vlsv")
   @compile_workload begin
      meta = load(initfile)
   end
end

end