module Vlasiator

using Requires
using StaticArrays: SVector, @SMatrix
using Printf: @sprintf
using LinearAlgebra: ×, dot, ⋅, norm, normalize!, normalize
using Statistics: mean
using EzXML
using Mmap: mmap
using WriteVTK
using LazyGrids: ndgrid
using LaTeXStrings
using Dates

include("utility/rotation.jl")
include("utility/log.jl")
include("utility/vector.jl")
include("utility/differential.jl")
include("vlsv/vlsvreader.jl")
include("vlsv/vlsvutility.jl")
include("utility/plot.jl")

export
   # vlsvreader
   MetaVLSV, VarInfo,
   load, readvariable, readparameter, readvariablemeta, readvcells,
   hasvariable, hasparameter, hasname, readmesh,
   # vlsvutility
   getcell, getslicecell, getlevel, refineslice, getcellcoordinates, getvcellcoordinates,
   getchildren, getparent, isparent, getsiblings,
   getcellinline, getnearestcellwithvdf, getcellwithvdf,
   getdensity, getvelocity, getpressure, getmaxwellianity, write_vtk, write_vlsv, issame,
   # plot helper
   SI, RE, Log, Linear, SymLog,
   # log
   readlog,
   # differential
   fg_grad, fg_curl, fg_div, fg_matder, fg_curvature, fg_normalize, fg_kappac, fg_kappa

precompile(load, (String,))
precompile(readvariable, (MetaVLSV, String, Bool))

# oo, a typedef!
realtype=Float64

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      include("plot/plots.jl")
   end
end

end
