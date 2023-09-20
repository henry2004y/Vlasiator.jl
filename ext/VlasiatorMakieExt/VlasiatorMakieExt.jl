module VlasiatorMakieExt

using Vlasiator, Printf
using Vlasiator: AxisUnit, ColorScale
using StatsBase: fit, Histogram, weights
using Makie.LaTeXStrings: latexstring

import Vlasiator: viz, viz!, vlheatmap, vlslice, vlslices, vdfvolume, vdfslice, vdfslices
import Makie

include("fullrecipe.jl")
include("typerecipe.jl")
include("interactive.jl")

end