var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"User Guide","title":"Examples","text":"","category":"section"},{"location":"examples/#Loading-VLSV-data","page":"User Guide","title":"Loading VLSV data","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Read meta data","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"filename = \"bulk.0000004.vlsv\"\nmeta = readmeta(filename)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Read variable meta data","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"showvariablemeta(meta, \"proton/vg_rho\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"A list of utility functions has been implemented for checking variable status. See here for the full list. ","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Read variable","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"data = readvariable(meta, \"proton/vg_rho\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"The same interface works for both DCCRG grid and FS grid variables. By default the returned DCCRG grid variable array is sorted by cell IDs. If in any case you want the original unsorted version as being stored in the file, simply say readvariable(meta, var, false).","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Get variable at a given location (This can be simplified even further later!)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"loc = [2.0, 0.0, 0.0]\nid = getcell(meta, loc)\nreadvariable(meta, \"proton/vg_rho\", id)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Get variable along a line between two points","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"point1 = [12Re, 0, 0]\npoint2 = [15Re, 0, 0]\ncellids, distances, coords = getcellinline(meta, point1, point2)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Compare VLSV files","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"One may want to check if two vlsv files are identical. This is tricky because","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"the structure of VLSV format does not guarantee parallel writing order;\nnumerical error accumulates with floating point representation.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"The key is that we should not check quantities that are related to MPI writing sequence. Note that even file sizes may vary depending on the number of MPI processes!","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"compare(filename1, filename2)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"There is an optional third argument to compare for setting the relative difference tolerance, with default being 1e-4. In practice relative difference works better for \"large\" numbers, and absolute difference works better for \"small\" numbers.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"More examples can be found in the repo.","category":"page"},{"location":"examples/#Computing-derived-quantities","page":"User Guide","title":"Computing derived quantities","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Vlasiator is capable of computing moments and some derived quantities and save them directly into VLSV files. More derived quantities computed from the saved quantities are also available in postprocessing, such as plasma β, velocity parallel/perpendicular to the magnetic field, pressure tensor with the third axis aligned with the magnetic field direction and so on. To avoid confusion about variable names, the convention here is that","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"if it is directly stored in the VLSV file, read the raw data;\notherwise check the availability in the derived variable list. All predefined names start with a capital letter.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To obtain a derived quantity, for example,","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"beta = readvariable(meta, \"Beta\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Here is a full list of available quantities[1]:","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Derived variable name Meaning Required variable[2]\nBmag magnetic field magnitude vg_b_vol\nEmag electric field magnitude vg_e_vol\nVmag bulk speed vg_v\nVS sound speed vg_ptensor_diagonal; vg_rho\nVA Alfvén speed vg_rho; Bmag\nMA Alfvén Mach number Vmag; VA\nVpar bulk velocity parallelmathbfB vg_v; vg_b_vol\nVperp bulk velocity perp mathbfB vg_v; vg_b_vol\nP scalar thermal pressure vg_ptensor_diagonal\nT scalar temperature P; vg_rho\nTpar temperature parallelmathbfB vg_rho; vg_ptensor_diagonal; vg_b_vol\nTperp temperature perp mathbfB vg_rho; vg_ptensor_offdiagonal; vg_b_vol\nProtated pressure tensor with widehatz parallel mathbfB vg_b_vol; vg_ptensor_diagonal; vg_ptensor_offdiagonal\nAnisotropy P_perp  P_parallel ptensor; B\nPdynamic dynamic pressure vg_rho; Vmag\nPoynting Poynting flux E; B\nBeta plasma beta, P  P_B P; vg_b_vol\nIonInertial ion inertial length vg_rho\nLarmor Larmor radius Vperp; Bmag\nGyrofrequency ion gyroperiod \nPlasmaperiod plasma oscillation period ","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"which can also be found as keys of dictionary in vlsvvariables.jl.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"[1]: For species specific variables, you need to add the species name at the front, separated by a slash. For example, the proton bulk velocity is a string proton/vg_v.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"[2]: If a required variable exists in the VLSV file, we try to use it directly instead of calculating from other variables. The interpolated FS grid variables onto DCCRG grid are preferred over original FS grid variables.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"warning: Warning\nThis part has not been carefully tested so it might not work or just generate wrong results!","category":"page"},{"location":"examples/#Plotting","page":"User Guide","title":"Plotting","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Vlasiator.jl does not include any plotting library as explicit dependency, but it offers plotting functionalities once the target plotting package is used.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Currently I would recommend using PyPlot.jl. Plots.jl is catching up, but it is still slower and lack of features. Makie.jl will be supported in the future if 3D plotting is necessary.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"More examples of customized plots can be found in the repo.","category":"page"},{"location":"examples/#PyPlot-Backend","page":"User Guide","title":"PyPlot Backend","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To trigger the Matplotlib plotting, using PyPlot. All the functions with identical names as in Matplotlib accept all possible keyword arguments supported by their Matplotlib counterparts, e.g. font width, font size, colormap, etc.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"note: Note\nThe method call to certain axes is not dispatched, e.g. ax.plot; as an alternative, one needs to pass ax as the third argument to the functions, e.g. plot(meta, \"rho\", ax)!","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Scalar colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"pcolormesh(meta, \"rho\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Vector z component colored contour for 2D simulation in a manually set range","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"pcolormesh(meta, \"rho\", op=:z, colorscale=Log, axisunit=RE, vmin=1e6, vmax=2e6)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Derived quantity colored contour for 2D simulation (as long as the input variable is in the predefined dictionary)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"pcolormesh(meta, \"b\", op=:z, colorscale=Linear, axisunit=SI)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Streamline for 2D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"streamplot(meta, \"rho_v\", comp=\"xy\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Quiver for 2D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"quiver(meta, \"rho_v\", comp=\"xy\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"The comp option is used to specify the two vector components.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"note: Note\nCurrently there is limited support for derived variables. This will be expanded and changed later for ease of use!","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"You can choose to use linear/log color scale via colorscale=Linear or colorscale=Log, plot vector components via e.g. op=:x or magnitude by default, and set unit via axisunit=RE etc..","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Mesh denoted by cell centers","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"plotmesh(meta; projection=\"z\", color=\"w\")","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Cut slice colored contour for 3D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"pcolormesh(meta, \"proton/vg_rho\", normal=:y, origin=0.0)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Velocity distribution function near a given spatial location coordinates = [0.0, 0.0, 0.0]","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"plot_vdf(meta, coordinates)","category":"page"},{"location":"examples/#Plots-Backend","page":"User Guide","title":"Plots Backend","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To trigger the Plots package plotting, using Plots. This backend supports all available attributes provided by Plots.jl. By default it uses GR, but a wide range of other options are also presented.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Scaler colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"heatmap(meta, \"rho\", aspect_ratio=:equal, c=:turbo)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Scaler colored contour with lines for 2D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"contourf(meta, \"rho)","category":"page"},{"location":"examples/#Gallery","page":"User Guide","title":"Gallery","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton density in advection flow with points denoting cell centers","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton density of Earth's magnetosphere in the meridional plane from 3D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton density of Earth's magnetosphere in the equatorial plane from 2D simulation, zoomed in to the magnetosheath and foreshock region, with streamlines and density contour at 1e7","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton density of Earth's magnetosphere in the meridional cut from 2D simulation, with fieldlines through fixed seeding points","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton density of Earth's magnetosphere in the normal cut planes from 3D simulation","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"Proton phase space distribution","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"(Image: )","category":"page"},{"location":"examples/#Converting-to-VTK","page":"User Guide","title":"Converting to VTK","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"We can convert VLSV files into VTK files! Since DCCRG is Cartesian based with uniform spacing, each level of refinement corresponds to a VTK image file, and the cell refinement relationships are defined by vtkGhostType as well as the vthb file.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To convert a VLSV file into VTK,","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"write_vtk(filename)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"This function accepts both file names and file meta.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To see the full list of options, please refer to the documentation in internal. Demo usage can be found here.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"warning: Warning\nAs of ParaView 5.9.1, there are display issues with VTKOverlappingAMR. However, we can read the generated image files directly. There is also an keyword argument for write_vtk called vti: when it is set to true, then only the image file at the highest refinement level is generated. This part is experimental and subject to change in the future.","category":"page"},{"location":"examples/#Calling-from-Python","page":"User Guide","title":"Calling from Python","text":"","category":"section"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"It is possible to call this package directly from Python with the aid of PyJulia. Following the installation steps described in the manual[3], and then inside Python REPL:","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"# Handling initialization issue for Conda\nfrom julia.api import Julia\njl = Julia(compiled_modules=False)\n\nfrom julia import Vlasiator\nfilename = \"bulk1.0001000.vlsv\"\nmeta = Vlasiator.readmeta(filename)\nvar = \"proton/vg_rho\"\ndata = Vlasiator.readvariable(meta, var)","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"To run a Julia script in Python,","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"# Handling initialization issue for Conda\nfrom julia.api import Julia\njl = Julia(compiled_modules=False)\njl.eval('include(\"examples/demo_2dplot_pyplot.jl\")')\nimport matplotlib.pyplot as plt\nplt.show()","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"note: Note\nThis approach is for you to have a taste of the package. For better integrated experience with its full power, I recommend using the package inside Julia.","category":"page"},{"location":"examples/","page":"User Guide","title":"User Guide","text":"[3]: For Debian-based Linux distributions, it gets a little bit tricky. Please refer to Troubleshooting for details.","category":"page"},{"location":"internal/#Internal","page":"API Reference","title":"Internal","text":"","category":"section"},{"location":"internal/#Public-APIs","page":"API Reference","title":"Public APIs","text":"","category":"section"},{"location":"internal/","page":"API Reference","title":"API Reference","text":"Modules = [Vlasiator]\nPrivate = false\nOrder = [:constant, :type, :function]","category":"page"},{"location":"internal/#Vlasiator.MetaData","page":"API Reference","title":"Vlasiator.MetaData","text":"Meta data declaration.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.VarInfo","page":"API Reference","title":"Vlasiator.VarInfo","text":"Variable metadata from the vlsv footer.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Base.ndims-Tuple{MetaData}","page":"API Reference","title":"Base.ndims","text":"ndims(meta) -> Int\n\nReturn the dimension of VLSV data.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.compare","page":"API Reference","title":"Vlasiator.compare","text":"compare(filename1, filename2, tol=1e-4) -> Bool\n\nCheck if two VLSV files are approximately identical.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.getcell-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.getcell","text":"getcell(meta, location) -> Int\n\nReturn cell ID containing the given spatial location.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getcellcoordinates-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.getcellcoordinates","text":"getcellcoordinates(meta, cellid) -> Vector{Float}\n\nReturn a given cell's coordinates.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getcellinline-Tuple{MetaData, Any, Any}","page":"API Reference","title":"Vlasiator.getcellinline","text":"getcellinline(meta, point1, point2) -> cellids, distances, coords\n\nReturns cell IDs, distances and coordinates for every cell in a line between two given points point1 and point2. May be improved later with preallocation!\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getchildren-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.getchildren","text":"getchildren(meta, cellid) -> Vector{Int}\n\nReturn direct children of cellid.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getlevel-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.getlevel","text":"getlevel(meta, cellid) -> Int\n\nReturn the AMR level of a given cell ID. Note that this function does not check if the VLSV file of meta actually contains cellid: it may be shadowed by refined children.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getnearestcellwithvdf-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.getnearestcellwithvdf","text":"getnearestcellwithvdf(meta, id) -> Int\n\nFind the nearest spatial cell with VDF saved of a given cell id in the file meta.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getparent-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.getparent","text":"getparent(meta, cellid) -> Int\n\nReturn the parent cell ID of given child cellid.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getsiblings-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.getsiblings","text":"getsiblings(meta, cellid) -> Vector{Int}\n\nReturn sibling cells of a given cellid, including itself.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getslicecell-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.getslicecell","text":"getslicecell(meta, slicelocation;\n   xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf) -> idlist, indexlist\n\nFind the cell ids idlist which are needed to plot a 2d cut through of a 3d mesh, in a direction with non infinity range at slicelocation, and the indexlist, which is a mapping from original order to the cut plane and can be used to select data onto the plane.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getvcellcoordinates-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.getvcellcoordinates","text":"getvcellcoordinates(meta, vcellids, pop=\"proton\")\n\nReturn velocity cells' coordinates of population pop and id vcellids.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasname-Tuple{Any, Any, Any}","page":"API Reference","title":"Vlasiator.hasname","text":"Check if the XMLElement elem contains a tag with name.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasparameter-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.hasparameter","text":"hasparameter(meta, param) -> Bool\n\nCheck if vlsv file contains a parameter.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasvariable-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.hasvariable","text":"hasvariable(meta, var) -> Bool\n\nCheck if the VLSV file contains a variable.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasvdf-Tuple{MetaData}","page":"API Reference","title":"Vlasiator.hasvdf","text":"hasvdf(meta) -> Bool\n\nCheck if VLSV file contains VDF.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.isparent-Tuple{MetaData, Integer}","page":"API Reference","title":"Vlasiator.isparent","text":"isparent(meta, cellid) -> Bool\n\nCheck if cellid is a parent cell.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readmeta-Tuple{AbstractString}","page":"API Reference","title":"Vlasiator.readmeta","text":"readmeta(filename; verbose=false) -> MetaData\n\nReturn MetaData from a vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readparameter-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.readparameter","text":"readparameter(meta, param)\n\nReturn the parameter value from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvariable","page":"API Reference","title":"Vlasiator.readvariable","text":"readvariable(meta, var, sorted=true) -> Array\n\nReturn variable value from the vlsv file. By default sorted=true, which means that for DCCRG grid the variables are sorted by cell ID.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.readvariable-Tuple{MetaData, AbstractString, Any}","page":"API Reference","title":"Vlasiator.readvariable","text":"readvariable(meta, var, ids) -> Array\n\nRead a variable var in a collection of cells ids.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvariablemeta-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.readvariablemeta","text":"readvariablemeta(meta, var) -> VarInfo\n\nReturn VarInfo about var in the vlsv file linked to meta.           \n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvcells-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.readvcells","text":"readvcells(meta, cellid; pop=\"proton\")\n\nRead velocity cells from a spatial cell of ID cellid, and return a map of velocity cell ids and corresponding value.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.refineslice-Tuple{MetaData, Any, Any, Any}","page":"API Reference","title":"Vlasiator.refineslice","text":"refineslice(meta, idlist, data, normal) -> Array\n\nGenerate scalar data on the finest refinement level given cellids idlist and variable data on the slice perpendicular to normal.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.write_vtk-Tuple{MetaData}","page":"API Reference","title":"Vlasiator.write_vtk","text":"write_vtk(meta::MetaData; kwargs...)\nwrite_vtk(filename; kwargs...)\n\nConvert VLSV file to VTK format.\n\nKeyword arguments\n\nvars=[\"\"]: select which variables to convert.\nascii=false: output stored in ASCII or compressed binary format.\nvti=false: generate image files on the highest refinement level only.\nverbose=false: display logs during conversion.\n\n\n\n\n\n","category":"method"},{"location":"internal/#PyPlot-helpers","page":"API Reference","title":"PyPlot helpers","text":"","category":"section"},{"location":"internal/","page":"API Reference","title":"API Reference","text":"Modules = [Vlasiator]\nPages   = [\"plot/pyplot.jl\"]","category":"page"},{"location":"internal/#Private-APIs","page":"API Reference","title":"Private APIs","text":"","category":"section"},{"location":"internal/","page":"API Reference","title":"API Reference","text":"Modules = [Vlasiator]\nPublic = false","category":"page"},{"location":"internal/#Vlasiator.AxisUnit","page":"API Reference","title":"Vlasiator.AxisUnit","text":"Axis unit type. Currently supported: SI, RE.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.ColorScale","page":"API Reference","title":"Vlasiator.ColorScale","text":"Color scales type for 2D plots. Currently supported: Log, Linear.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.MeshInfo","page":"API Reference","title":"Vlasiator.MeshInfo","text":"Mesh size information.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.fillmesh-Tuple{MetaData, Any}","page":"API Reference","title":"Vlasiator.fillmesh","text":"fillmesh(meta::MetaData, vars; verbose=false)\n\nFill the DCCRG mesh with quantity of vars on all refinement levels.\n\nReturn arguments\n\ncelldata::Vector{Vector{Array}}: data for each variable on each AMR level.\nvtkGhostType::Array{UInt8}: cell status (to be completed!). \n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get1stcell-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.get1stcell","text":"Return the first cellid - 1 on my level.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getObjInfo-NTuple{5, Any}","page":"API Reference","title":"Vlasiator.getObjInfo","text":"Return size and type information for the object.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getRotationB-Tuple{Any}","page":"API Reference","title":"Vlasiator.getRotationB","text":"getRotationB(B) -> SMatrix\n\nObtain a rotation matrix with each column being a unit vector which is parallel (v3) and perpendicular (v1,v2) to the magnetic field B. The two perpendicular directions are chosen based on the reference vector of z-axis in the Cartesian coordinates.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_rotation_matrix-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.get_rotation_matrix","text":"get_rotation_matrix(vector, angle)\n\nCreates a rotation matrix that rotates around a unit vector by an angle in radians. References: https://en.wikipedia.org/wiki/Rodrigues'rotationformula https://en.wikipedia.org/wiki/Rotationmatrix#Rotationmatrixfromaxisandangle\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getfooter-Tuple{Any}","page":"API Reference","title":"Vlasiator.getfooter","text":"Return the xml footer of vlsv.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getindexes-NTuple{5, Any}","page":"API Reference","title":"Vlasiator.getindexes","text":"Compute every cell id's x, y and z indexes on the given refinement level (0-based).\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readmesh-NTuple{4, Any}","page":"API Reference","title":"Vlasiator.readmesh","text":"Return mesh related variable.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvector-NTuple{4, Any}","page":"API Reference","title":"Vlasiator.readvector","text":"Return vector data from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateTensorToVectorZ!-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.rotateTensorToVectorZ!","text":"rotateTensorToVector(tensor, vector)\n\nRotates tensor with a rotation matrix that aligns z-axis with vector.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateWithB!-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.rotateWithB!","text":"rotateWithB!(T, B)\n\nRotate the tensor T with the 3rd direction aligned with B. See also: rotateWithB\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateWithB-Tuple{Any, Any}","page":"API Reference","title":"Vlasiator.rotateWithB","text":"rotateWithB(T, B) -> Matrix\n\nRotate the tensor T with the 3rd direction aligned with B. See also: rotateWithB!\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.save_image","page":"API Reference","title":"Vlasiator.save_image","text":"save_image(meta::MetaData, file, vars, data, vtkGhostType, level,\n   xcells, ycells, zcells, ascii=false, append=true)\n\nSave data of name vars at AMR level into VTK image file of name file.\n\nArguments\n\nfile::String: output file name.\nvars::Vector{String}: variable names to be saved.\ndata::Vector{Vector}: data for all the variables on each refinement level.\nvtkGhostType::Array{UInt8}: array for visibility control.\nlevel::Int: refinement level (0-based).\nxcells, ycells, zcells: original mesh sizes.\nascii=false: save output in ASCII or binary format.\nappend=true: determines whether to append data at the end of file or do in-block writing.\n\n\n\n\n\n","category":"function"},{"location":"log/#Log","page":"Log","title":"Log","text":"","category":"section"},{"location":"log/#Performance","page":"Log","title":"Performance","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"The VLSV loader inherits the basic structure from Analysator and is redesigned for performance.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Besides the language difference in speed, one of the key decisions in boosting performance is to avoid the usage of dictionary with integer keys as much as possible.\nIt is generally faster to read a bunch of cell IDs together than to read each cell one-by-one.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"For development, it is recommended to use PkgBenchmark.jl to run the test suite:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"using PkgBenchmark, Vlasiator\nbenchmarkpkg(Vlasiator)","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"or if you want to compare the current status of the package against a different git version","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"judge(Vlasiator, \"97e3dca6b2474d7bdc5b62b5bf98ecf070516e5e\")","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"See more in the PkgBenchmark manual.","category":"page"},{"location":"log/#Benchmarks","page":"Log","title":"Benchmarks","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"Initial tests on reading variables from sample VLSV files: ","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"DCCRG grid","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Julia tmean [μs] Python tmean [μs]\n2MB 200 2MB 1000\n50MB 400 50MB 1000","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Field solver grid[1]","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 13\nPython 45","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"[1]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read fg_b once; it will go out of memory for the second time! This reading time corresponds to 35% of the maximum sequential read speed on the target machine.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Timing from starting Julia/Python to the first plot[2]: | 2.3GB  | tmean [s] | |:–––-|:––––-:| | Julia  | 11.6  | | Python | 9.3   |","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"[2]: This inefficieny is a famous problem in Julia known as \"time to first plot\". On the Python side, however, I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s).","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Reading and plotting one 2d slice of proton density out of 3D AMR data:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 0.35\nPython 1.7","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Virtual satellite tracking from 845 frames of 3D AMR data (26G per frame) on Vorna:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"1 CPU tmean [m][3]\nJulia 11\nPython 125","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"[3]: The timings are for a single CPU. With multithreading, the Julia timings can scale linearly on a node with the number of cores used. For example, with 8 threads, Julia takes ~80s to finish.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Vlasiator","category":"page"},{"location":"#Vlasiator.jl","page":"Home","title":"Vlasiator.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics Vlasiator. This lightweight package is built upon its sister in Python Analysator, and is carefully designed for performance and ease of use. It becomes even more powerful when using together with other fantastic packages: combined with external packages like FieldTracer.jl and TestParticle.jl, it is possible to do all kinds of in-depth analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Vlasiator.jl contains the following features:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Reading VLSV format data.\nConverting VLSV into VTK format.\nExtracting quantities from the simulation at a given point/line/cut.\nPlotting 2D cuts.\nPlotting phase space distributions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThis package is still young, so be careful for any future breaking changes!","category":"page"},{"location":"","page":"Home","title":"Home","text":"warning: Warning\nThis package mostly aims at supporting Vlasiator 5.0+. Older versions of Vlasiator has different naming standard for outputs, and is not guaranteed to work. This analysator wiki page describes the old and new naming standards in detail.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install it,","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can then get started with","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to use Plots.jl for visualization, add it also through the pkg manager; if you aim at using Matplotlib, besides adding PyPlot, you should also link to a preinstalled Python version by setting the environment variable and building the PyCall package","category":"page"},{"location":"","page":"Home","title":"Home","text":"ENV[\"PYTHON\"]=\"your python executable\"\nPkg.build(\"PyCall\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Details are described in automated matplotlib installation.","category":"page"},{"location":"#Author","page":"Home","title":"Author","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This module is written by Hongyang Zhou.","category":"page"}]
}
