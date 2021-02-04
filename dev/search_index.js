var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Example","title":"Examples","text":"","category":"section"},{"location":"examples/#Loading-VLSV-data","page":"Example","title":"Loading VLSV data","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Read meta data","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"filename = \"bulk.0000004.vlsv\"\nmeta = read_meta(filename)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Display all variable names","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"vars = show_variables(meta)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Read variable","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"data = read_variable(meta, \"proton/vg_rho\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"By default the returned variable array is sorted by cell IDs. If in any case you want the original unsorted version as being stored in the file, simply say read_variable(meta, var, false).","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Get variable at a given location (This can be simplified even further later!)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"loc = [2.0, 0.0, 0.0]\nid = get_cellid(meta, loc)\nread_variable_select(meta, \"proton/vg_rho\", id)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Get variable along a line between two points","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"point1 = [12Re, 0, 0]\npoint2 = [15Re, 0, 0]\ncellids, distances, coords = get_cell_in_line(meta, point1, point2)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Combined with external packages like FieldTracer.jl, it is possible to do all kinds of in-depth analysis. More examples can be found in the repo.","category":"page"},{"location":"examples/#Plotting","page":"Example","title":"Plotting","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Vlasiator.jl does not have any plotting library as dependency, but it offers plotting functionalities through additional scripts. To use a specific plotting library, just include the target script (e.g. pyplot.jl) under src/plot folder:","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"include(\"src/plot/pyplot.jl\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Currently I would recommend using PyPlot.jl. Plots.jl is catching up, but it is still slower and lack of features. Makie.jl will be supported in the future if 3D plotting is necessary.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"More examples of customized plots can be found in the repo.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Sample outputs:","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the meridional plane from 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the equatorial plane from 2D simulation, zoomed in to the magnetosheath and foreshock region, with streamlines and density contour at 1e7","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the normal cut planes from 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/#PyPlot-Backend","page":"Example","title":"PyPlot Backend","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Scaler colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plot_pcolormesh(meta, \"rho\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Vector z component colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plot_pcolormesh(meta, \"rho\", op=\"z\", islinear=false, axisunit=\"Re\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Derived quantity colored contour for 2D simulation (as long as the input variable is in the predefined dictionary)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plot_pcolormesh(meta, \"b\", op=\"z\", islinear=false, axisunit=\"Re\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Streamline for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"streamline(meta, \"rho_v\", comp=\"xy\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"The comp option is used to specify the two vector components.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"note: Note\nCurrently there is limited support for derived variables. This will be expanded and changed later for ease of use!","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"You can choose to use linear/log color scale, plot vector components via e.g. op=\"x\" or magnitude by default, and set unit via axisunit=\"Re\" etc..","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Cut slice colored contour for 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plot_colormap3dslice(meta, \"proton/vg_rho\", normal=\"y\")","category":"page"},{"location":"examples/#Plots.jl-Backend","page":"Example","title":"Plots.jl Backend","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Scaler colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"heatmap(meta, \"rho\", aspect_ratio=:equal, c=:turbo)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Scaler colored contour with lines for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"contourf(meta, \"rho)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"This backend supports all available attributes provided by Plots.jl.","category":"page"},{"location":"internal/#Internal","page":"Internal","title":"Internal","text":"","category":"section"},{"location":"internal/#Public-APIs","page":"Internal","title":"Public APIs","text":"","category":"section"},{"location":"internal/","page":"Internal","title":"Internal","text":"Modules = [Vlasiator]\nPrivate = false\nOrder = [:function]","category":"page"},{"location":"internal/#Vlasiator.getSliceCellID-Tuple{Any,Any,Any}","page":"Internal","title":"Vlasiator.getSliceCellID","text":"getSliceCellID(meta, slicelocation, maxreflevel; xmin=-Inf, xmax=Inf,\n   ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf)\n\nFind the cell ids which are needed to plot a 2d cut through of a 3d mesh.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_amr_level-Tuple{Any,Any}","page":"Internal","title":"Vlasiator.get_amr_level","text":"get_amr_level(meta::MetaData, cellid)\n\nReturn the AMR level of a given cell id.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_cell_coordinates-Tuple{Any,Any}","page":"Internal","title":"Vlasiator.get_cell_coordinates","text":"get_cell_coordinates(meta::MetaData, cellid)\n\nReturn a given cellid's coordinates.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_cell_in_line-Tuple{Any,Any,Any}","page":"Internal","title":"Vlasiator.get_cell_in_line","text":"get_cell_in_line(meta::MetaData, point1, point2)\n\nReturns cell IDs, distances and coordinates for every cell in a line between two given points. May be improved later with preallocation!\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_cellid-Tuple{Any,Any}","page":"Internal","title":"Vlasiator.get_cellid","text":"get_cell_id(meta::MetaData, location)\n\nReturn cell ID containing the given spatial location.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_max_amr_level-Tuple{Any}","page":"Internal","title":"Vlasiator.get_max_amr_level","text":"get_max_amr_level(meta::MetaData)\n\nFind the highest refinement level.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.has_name-Tuple{Any,Any,Any}","page":"Internal","title":"Vlasiator.has_name","text":"Check if the XMLElement elem contains a tag with name.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.has_parameter-Tuple{MetaData,Any}","page":"Internal","title":"Vlasiator.has_parameter","text":"has_parameter(meta::MetaData, param)\n\nCheck if vlsv file contains a parameter.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.has_variable-Tuple{Any,Any}","page":"Internal","title":"Vlasiator.has_variable","text":"Check if the VLSV file contains a variable.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.read_meta-Tuple{AbstractString}","page":"Internal","title":"Vlasiator.read_meta","text":"read_meta(filename; verbose=false)\n\nReturn a struct of MetaData from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.read_parameter-Tuple{MetaData,Any}","page":"Internal","title":"Vlasiator.read_parameter","text":"read_parameter(meta::MetaData, param)\n\nReturn the parameter value from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.read_variable","page":"Internal","title":"Vlasiator.read_variable","text":"read_variable(meta::MetaData, var, sorted=true)\n\nReturn variable value from the vlsv file. For DCCRG grid, the variables are sorted by cell ID by default.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.read_variable_info-Tuple{Any,Any}","page":"Internal","title":"Vlasiator.read_variable_info","text":"read_variable_info(meta::MetaData, var)\n\nReturn a struct of VarInfo.           \n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.read_variable_select","page":"Internal","title":"Vlasiator.read_variable_select","text":"read_variable_select(meta::MetaData, var, idIn=UInt[])\n\nRead a variable in a collection of cells.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.refine_data-NTuple{5,Any}","page":"Internal","title":"Vlasiator.refine_data","text":"refine_data(meta, idlist, data, maxreflevel, normal)\n\nGenerate scalar data on the finest refinement level.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.show_variables-Tuple{MetaData}","page":"Internal","title":"Vlasiator.show_variables","text":"Display all variables in the VLSV file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Private-types-and-methods","page":"Internal","title":"Private types and methods","text":"","category":"section"},{"location":"internal/","page":"Internal","title":"Internal","text":"Modules = [Foo]\nPublic = false","category":"page"},{"location":"log/#Log","page":"Log","title":"Log","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"This package was born when I was learning Vlasiator and its corresponding data structures. The VLSV loader inherits the basic structure from Analysator and is redesigned for performance.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"The function APIs are made to be consistent with Analysator.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"The IOstream handle for the VLSV file requires some special attention. In the current implementation, once the meta data is read, the file stays open until one explictly says close(meta.fid). On the Windows platform, it is not allowed to delete the file before the IO is closed. However, this is allowed in Unix, so be careful.","category":"page"},{"location":"log/#Plotting-Philosophy","page":"Log","title":"Plotting Philosophy","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"We should not take over what underlying plotting libraries like Matplotlib offers. Users should be able to modify the figures as they wish even if they only know how to use the well-known plotting libraries. Do not reinvent the wheel. For customized plotting, simply provide some sample scripts for the common tasks like zooming-in, change font sizes, add text boxes, etc..","category":"page"},{"location":"log/#Benchmarks","page":"Log","title":"Benchmarks","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"Initial tests on reading variables from a VLSV file: ","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"DCCRG grid","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"2MB tmean [μs]\nJulia 200\nPython 1000","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"50MB tmean [μs]\nJulia 400\nPython 1000","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Field solver grid[1]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant. With 16 GB memory it is barely enough to read fg_b once; it will go out of memory for the second time!","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 13\nPython 45","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s). Same thing for Julia costs 1.0s (first time ~8s including everything).","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Reading and plotting one 2d slice of proton density out of 3D AMR data:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 0.35\nPython 1.7","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Virtual satellite tracking from 845 frames of 3D AMR data (26G per frame) on Vorna:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"1 CPU tmean [m]\nJulia 11\nPython 125","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Note that the above timings are for a single CPU. With only one command added for multithreading, the Julia timings can be improved by n where n is the number of threads. For example, with 8 threads, Julia takes ~80s to finish.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Vlasiator","category":"page"},{"location":"#Vlasiator","page":"Home","title":"Vlasiator","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics Vlasiator. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package contains the following features:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Reading VLSV format data.\nExtracting quantities from the simulation at a given point/line/cut.\nPlotting 2D cuts from 2D/3D simulations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThis package is still young, so be careful for any future breaking changes!","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install it,","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can then get started with","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"More usages can be found in the examples.","category":"page"}]
}
