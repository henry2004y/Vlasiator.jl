# Sample script for plotting the phase space density near a given spatial
# location.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot

filename = "bulk.0000001.vlsv"

meta = readmeta(filename)

coordinates = [0.0, 0.0, 0.0]

plot_vdf(meta, coordinates; verbose=true)