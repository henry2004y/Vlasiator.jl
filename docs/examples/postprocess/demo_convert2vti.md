# ---
# title: Converting to VTK
# id: demo_convertvtk
# date: 2023-02-25
# author: "[Hongyang Zhou](https://github.com/henry2004y)"
# julia: 1.8.5
# description: This demo shows how to convert VLSV to VTK
# ---

To convert VLSV time series files to VTK image files using multithreads,
```julia
using Vlasiator, Glob

files = glob("*.vlsv")

Threads.@threads for file in files
   @info file, Threads.threadid()
   write_vtk(file; vti=true)
end
```