# Sample script for converting VLSV time series files to VTK image files.
# For multithreading run,
#   julia -t 2 demo_convert2vti.jl
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, Glob

files = glob("*.vlsv")

Threads.@threads for i in 1:length(files)
   file = files[i]
   @info file, Threads.threadid()
   write_vtk(file; vti=true)
end