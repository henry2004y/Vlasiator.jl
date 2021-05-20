# Sample script for converting VLSV time series files to VTK image files.
# For multithreading run,
#   julia -t 2 demo_convert2vti.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob

filenames = glob("*.vlsv")

Threads.@threads for i in 1:length(filenames)
   fname = filenames[i]
   @info fname, Threads.threadid()
   write_vtk(fname; vti=true)
end