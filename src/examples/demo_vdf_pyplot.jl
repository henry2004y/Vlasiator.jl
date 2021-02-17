# Sample script for plotting the phase space density near a given spatial
# location.
#
# Hongyang Zhou, hyzhou@umich.edu 02/17/2021

include("../plot/pyplot.jl")

filename = "bulk.0000004.vlsv"

meta = read_meta(filename)

coordinates = [0.0, 0.0, 0.0]

plot_vdf(meta, coordinates; verbose=true)