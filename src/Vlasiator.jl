module Vlasiator

# Hongyang Zhou, hyzhou@umich.edu

using Requires

include("vlsv/vlsvreader.jl")
include("vlsv/vlsvutility.jl")

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("plot/pyplot.jl")
   @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plot/plots.jl")
end

end
