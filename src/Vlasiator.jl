module Vlasiator

# Hongyang Zhou, hyzhou@umich.edu

using Requires

include("vlsv/vlsvreader.jl")
include("vlsv/vlsvutility.jl")
include("plot/plots.jl")

function __init__()
   @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("plot/pyplot.jl")
end

end
