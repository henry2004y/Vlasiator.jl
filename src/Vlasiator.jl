module Vlasiator

using Requires, UnPack

include("utility/rotation.jl")
include("utility/plot.jl")
include("utility/log.jl")
include("utility/curvature.jl")
include("vlsv/vlsvreader.jl")
include("vlsv/vlsvutility.jl")

precompile(load, (String,))

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
      include("plot/pyplot.jl")
   end
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      include("plot/plots.jl")
   end
end

end
