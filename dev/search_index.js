var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Example","title":"Examples","text":"","category":"section"},{"location":"examples/#Loading-VLSV-data","page":"Example","title":"Loading VLSV data","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Read meta data","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"filename = \"bulk.0000004.vlsv\"\nmeta = readmeta(filename)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Display all variable names","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"vars = showvariables(meta)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Read variable meta data","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"showvariablemeta(meta, \"proton/vg_rho\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"A list of utility functions has been implemented for checking variable status. See here for the full list. ","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Read variable","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"data = readvariable(meta, \"proton/vg_rho\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"The same interface works for both DCCRG grid and FS grid variables. By default the returned DCCRG grid variable array is sorted by cell IDs. If in any case you want the original unsorted version as being stored in the file, simply say readvariable(meta, var, false).","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Get variable at a given location (This can be simplified even further later!)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"loc = [2.0, 0.0, 0.0]\nid = getcell(meta, loc)\nreadvariable(meta, \"proton/vg_rho\", id)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Get variable along a line between two points","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"point1 = [12Re, 0, 0]\npoint2 = [15Re, 0, 0]\ncellids, distances, coords = getcellinline(meta, point1, point2)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Compare VLSV files","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"One may want to check if two vlsv files are identical. This is tricky because","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"the structure of VLSV format does not guarantee parallel writing order;\nnumerical error accumulates with floating point representation.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"The key is that we should not check quantities that are related to MPI writing sequence. Note that even file sizes may vary depending on the number of MPI processes!","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"compare(filename1, filename2)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"There is an optional third argument to compare for setting the relative difference tolerance, with default being 1e-4. In practice relative difference works better for \"large\" numbers, and absolute difference works better for \"small\" numbers.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"More examples can be found in the repo.","category":"page"},{"location":"examples/#Computing-derived-quantities","page":"Example","title":"Computing derived quantities","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Vlasiator is capable of computing moments and some derived quantities and save them directly into VLSV files. More derived quantities computed from the saved quantities are also available in postprocessing, such as plasma β, velocity parallel/perpendicular to the magnetic field, pressure tensor with the third axis aligned with the magnetic field direction and so on. To avoid confusion about variable names, the convention here is that","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"if it is directly stored in the VLSV file, read the raw data;\notherwise check the availability in the derived variable list. All predefined names start with a capital letter.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"To obtain a derived quantity, for example,","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"beta = readvariable(meta, \"Beta\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Here is a full list of available quantities[1]:","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Derived variable name Meaning Required variable[2]\nBmag magnetic field magnitude vg_b_vol\nEmag electric field magnitude vg_e_vol\nVmag bulk speed vg_v\nVS sound speed vg_ptensor_diagonal; vg_rho\nVA Alfvén speed vg_rho; Bmag\nMA Alfvén Mach number Vmag; VA\nVpar bulk velocity parallelmathbfB vg_v; vg_b_vol\nVperp bulk velocity perp mathbfB vg_v; vg_b_vol\nP scalar thermal pressure vg_ptensor_diagonal\nT scalar temperature P; vg_rho\nTpar temperature parallelmathbfB vg_rho; vg_ptensor_diagonal; vg_b_vol\nTperp temperature perp mathbfB vg_rho; vg_ptensor_offdiagonal; vg_b_vol\nProtated pressure tensor with widehatz parallel mathbfB vg_b_vol; vg_ptensor_diagonal; vg_ptensor_offdiagonal\nAnisotropy P_perp  P_parallel ptensor; B\nPdynamic dynamic pressure vg_rho; Vmag\nPoynting Poynting flux E; B\nBeta plasma beta, P  P_B P; vg_b_vol\nIonInertial ion inertial length vg_rho\nLarmor Larmor radius Vperp; Bmag\nGyrofrequency ion gyroperiod \nPlasmaperiod plasma oscillation period ","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"which can also be found as keys of dictionary in vlsvvariables.jl.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"[1]: For species specific variables, you need to add the species name at the front, separated by a slash. For example, the proton bulk velocity is a string proton/vg_v.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"[2]: If a required variable exists in the VLSV file, we try to use it directly instead of calculating from other variables. The interpolated FS grid variables onto DCCRG grid are preferred over original FS grid variables.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"warning: Warning\nThis part has not been carefully tested so it might not work or just generate wrong results!","category":"page"},{"location":"examples/#Plotting","page":"Example","title":"Plotting","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Vlasiator.jl does not include any plotting library as explicit dependency, but it offers plotting functionalities once the target plotting package is used.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Currently I would recommend using PyPlot.jl. Plots.jl is catching up, but it is still slower and lack of features. Makie.jl will be supported in the future if 3D plotting is necessary.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"More examples of customized plots can be found in the repo.","category":"page"},{"location":"examples/#PyPlot-Backend","page":"Example","title":"PyPlot Backend","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"To trigger the Matplotlib plotting, using PyPlot. All the functions with identical names as in Matplotlib accept all possible keyword arguments.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"note: Note\nThe method call to certain axes is not dispatched, e.g. ax.plot; as an alternative, one needs to pass ax as the third argument to the functions, e.g. plot(meta, \"rho\", ax)!","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Scalar colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"pcolormesh(meta, \"rho\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Vector z component colored contour for 2D simulation in a manually set range","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"pcolormesh(meta, \"rho\", op=:z, colorscale=Log, axisunit=RE, vmin=1e6, vmax=2e6)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Derived quantity colored contour for 2D simulation (as long as the input variable is in the predefined dictionary)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"pcolormesh(meta, \"b\", op=:z, colorscale=Linear, axisunit=SI)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Streamline for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"streamplot(meta, \"rho_v\", comp=\"xy\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Quiver for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"quiver(meta, \"rho_v\", comp=\"xy\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"The comp option is used to specify the two vector components.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"note: Note\nCurrently there is limited support for derived variables. This will be expanded and changed later for ease of use!","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"You can choose to use linear/log color scale via colorscale=Linear or colorscale=Log, plot vector components via e.g. op=:x or magnitude by default, and set unit via axisunit=RE etc..","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Mesh denoted by cell centers","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plotmesh(meta; projection=\"z\", color=\"w\")","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Cut slice colored contour for 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"pcolormesh(meta, \"proton/vg_rho\", normal=:y, origin=0.0)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Velocity distribution function near a given spatial location coordinates = [0.0, 0.0, 0.0]","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"plot_vdf(meta, coordinates)","category":"page"},{"location":"examples/#Plots-Backend","page":"Example","title":"Plots Backend","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"To trigger the Plots package plotting, using Plots. This backend supports all available attributes provided by Plots.jl. By default it uses GR, but a wide range of other options are also presented.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Scaler colored contour for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"heatmap(meta, \"rho\", aspect_ratio=:equal, c=:turbo)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Scaler colored contour with lines for 2D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"contourf(meta, \"rho)","category":"page"},{"location":"examples/#Gallery","page":"Example","title":"Gallery","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density in advection flow with points denoting cell centers","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the meridional plane from 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the equatorial plane from 2D simulation, zoomed in to the magnetosheath and foreshock region, with streamlines and density contour at 1e7","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the meridional cut from 2D simulation, with fieldlines through fixed seeding points","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton density of Earth's magnetosphere in the normal cut planes from 3D simulation","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"Proton phase space distribution","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"(Image: )","category":"page"},{"location":"examples/#Calling-from-Python","page":"Example","title":"Calling from Python","text":"","category":"section"},{"location":"examples/","page":"Example","title":"Example","text":"It is possible to call this package directly from Python with the aid of PyJulia. Following the installation steps described in the manual[2], and then inside Python REPL:","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"# Handling initialization issue for Conda\nfrom julia.api import Julia\njl = Julia(compiled_modules=False)\n\nfrom julia import Vlasiator\nfilename = \"bulk1.0001000.vlsv\"\nmeta = Vlasiator.readmeta(filename)\nvar = \"proton/vg_rho\"\ndata = Vlasiator.readvariable(meta, var)","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"To run a Julia script in Python,","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"# Handling initialization issue for Conda\nfrom julia.api import Julia\njl = Julia(compiled_modules=False)\njl.eval('include(\"examples/demo_2dplot_pyplot.jl\")')\nimport matplotlib.pyplot as plt\nplt.show()","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"note: Note\nThis approach is for you to have a taste of the package. For better integrated experience with its full power, I recommend using the package inside Julia.","category":"page"},{"location":"examples/","page":"Example","title":"Example","text":"[2]: For Debian-based Linux distributions, it gets a little bit tricky. Please refer to Troubleshooting for details.","category":"page"},{"location":"internal/#Internal","page":"Internal","title":"Internal","text":"","category":"section"},{"location":"internal/#Public-APIs","page":"Internal","title":"Public APIs","text":"","category":"section"},{"location":"internal/","page":"Internal","title":"Internal","text":"Modules = [Vlasiator]\nPrivate = false\nOrder = [:constant, :type, :function]","category":"page"},{"location":"internal/#Vlasiator.MetaData","page":"Internal","title":"Vlasiator.MetaData","text":"Meta data declaration.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.VarInfo","page":"Internal","title":"Vlasiator.VarInfo","text":"Variable metadata from the vlsv footer.\n\n\n\n\n\n","category":"type"},{"location":"internal/#PyPlot.pcolormesh","page":"Internal","title":"PyPlot.pcolormesh","text":"pcolormesh(meta::MetaData, var, ax=nothing;\n   op=:mag, axisunit=RE, colorscale=Log, vmin=-Inf, vmax=Inf, addcolorbar=true)\n\nPlot a variable using pseudocolor from 2D VLSV data. If ax is provided, then it will plot on that axes. If 3D or AMR grid detected, it will pass arguments to pcolormeshslice.\n\nOptional arguments\n\nop::Symbol: the component of a vector to plot, chosen from :mag, :x, :y, :z.\naxisunit::AxisUnit: the unit of axis ∈ RE, SI.\ncolorscale::ColorScale: whether to use linear scale for data.\nvmin::Float: minimum data range. Set to maximum of data if not specified. \nvmax::Float: maximum data range. Set to minimum of data if not specified.\naddcolorbar::Bool: whether to add a colorbar to the colormesh.\n\npcolormesh(meta, var)\n\npcolormesh(meta, var, axisunit=SI)\n\npcolormesh(data, func, colorscale=Linear)\n\n\n\n\n\n","category":"function"},{"location":"internal/#PyPlot.plot","page":"Internal","title":"PyPlot.plot","text":"plot(meta, var, ax=nothing; kwargs)\n\nPlot var from meta of 1D VLSV data. If ax===nothing, plot on the current active axes.\n\n\n\n\n\n","category":"function"},{"location":"internal/#PyPlot.quiver","page":"Internal","title":"PyPlot.quiver","text":"quiver(meta, var, ax=nothing; comp=\"xy\", axisunit=RE, kwargs...)\n\nWrapper over Matplotlib's quiver function. If ax===nothing, plot on the current active axes. The comp option can take a subset of \"xyz\" in any order. axisunit can be chosen from RE, SI. The keyword arguments can be any valid Matplotlib arguments into quiver.\n\n\n\n\n\n","category":"function"},{"location":"internal/#PyPlot.streamplot","page":"Internal","title":"PyPlot.streamplot","text":"streamplot(meta, var, ax=nothing; comp=\"xy\", axisunit=RE, kwargs...)\n\nWrapper over Matplotlib's streamplot function. The comp option can take a subset of \"xyz\" in any order. axisunit can be chosen from RE, SI. The keyword arguments can be any valid Matplotlib arguments into streamplot.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.compare","page":"Internal","title":"Vlasiator.compare","text":"compare(filename1, filename2, tol=1e-4) -> Bool\n\nCheck if two VLSV files are approximately identical.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.getcell-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.getcell","text":"getcell(meta, location) -> Int\n\nReturn cell ID containing the given spatial location.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getcellcoordinates-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.getcellcoordinates","text":"getcellcoordinates(meta, cellid) -> Vector{Float}\n\nReturn a given cell's coordinates.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getcellinline-Tuple{Any, Any, Any}","page":"Internal","title":"Vlasiator.getcellinline","text":"getcellinline(meta, point1, point2) -> cellids, distances, coords\n\nReturns cell IDs, distances and coordinates for every cell in a line between two given points point1 and point2. May be improved later with preallocation!\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getlevel-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.getlevel","text":"getlevel(meta, cellid) -> Int\n\nReturn the AMR level of a given cell ID.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getmaxamr-Tuple{Any}","page":"Internal","title":"Vlasiator.getmaxamr","text":"getmaxamr(meta) -> Int\n\nFind the highest refinement level of a given vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getnearestcellwithvdf-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.getnearestcellwithvdf","text":"getnearestcellwithvdf(meta, id) -> Int\n\nFind the nearest spatial cell with VDF saved of a given cell id in the file meta.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getslicecell-Tuple{Any, Any, Any}","page":"Internal","title":"Vlasiator.getslicecell","text":"getslicecell(meta, slicelocation, maxreflevel;\n   xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf) -> idlist, indexlist\n\nFind the cell ids idlist which are needed to plot a 2d cut through of a 3d mesh, in a direction with non infinity range at slicelocation, and the indexlist, which is a mapping from original order to the cut plane and can be used to select data onto the plane.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getvcellcoordinates-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.getvcellcoordinates","text":"getvcellcoordinates(meta, vcellids, pop=\"proton\")\n\nReturn velocity cells' coordinates of population pop and id vcellids.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasname-Tuple{Any, Any, Any}","page":"Internal","title":"Vlasiator.hasname","text":"Check if the XMLElement elem contains a tag with name.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasparameter-Tuple{MetaData, Any}","page":"Internal","title":"Vlasiator.hasparameter","text":"hasparameter(meta, param) -> Bool\n\nCheck if vlsv file contains a parameter.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasvariable-Tuple{MetaData, Any}","page":"Internal","title":"Vlasiator.hasvariable","text":"hasvariable(meta, var) -> Bool\n\nCheck if the VLSV file contains a variable.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.hasvdf-Tuple{MetaData}","page":"Internal","title":"Vlasiator.hasvdf","text":"hasvdf(meta) -> Bool\n\nCheck if VLSV file contains VDF.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.pcolormeshslice","page":"Internal","title":"Vlasiator.pcolormeshslice","text":"pcolormeshslice(meta, var, ax=nothing; kwargs...)\n\nPlot pseudocolor var on a 2D slice of 3D vlsv data. If ax is provided, then it will plot on that axes. It would be easier to call pcolormesh.\n\nOptional arguments\n\nop::Symbol: the component of a vector to plot, chosen from :mag, :x, :y, :z.\norigin::Float: center of slice plane in the normal direction.\nnormal::Symbol: the normal direction of cut plane, chosen from :x, :y, :z.\naxisunit::AxisUnit: the unit of axis ∈ RE, SI.\ncolorscale::ColorScale: color scale for data ∈ (Linear, Log)\nvmin::Real: minimum data range. Set to maximum of data if not specified. \nvmax::Real: maximum data range. Set to minimum of data if not specified.\naddcolorbar::Bool: whether to add a colorbar to the colormesh.\n\npcolormeshslice(meta, var)\n\npcolormeshslice(meta, var, op=:z, origin=1.0, normal=:x)\n\npcolormeshslice(data, func, colorscale=Log)\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.plot_vdf","page":"Internal","title":"Vlasiator.plot_vdf","text":"plot_vdf(meta, location, ax=nothing; kwargs...)\n\nPlot the 2D slice cut of phase space distribution function at location within velocity range limits. If ax===nothing, plot on the current active axes.\n\nOptional arguments\n\nunit::AxisUnit: axis unit in SI, RE.\nlimits::Vector{Real}: velocity space range given in [xmin, xmax, ymin, ymax].\nslicetype: string for choosing the slice type from \"xy\", \"xz\", \"yz\", \"bperp\", \"bpar\",\n\n\"bpar1\".\n\ncenter: string for setting the reference frame from \"bulk\", \"peak\".\nvslicethick: setting the velocity space slice thickness in the normal direction. If set\n\nto 0, the whole distribution along the normal direction is projected onto a plane. Currently this is only meaningful when center is set such that a range near the bulk/peak normal velocity is selected! \n\nweight::Symbol: choosing distribution weights from phase space density or particle flux\n\nbetween :particle and :flux.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.plotmesh","page":"Internal","title":"Vlasiator.plotmesh","text":"plotmesh(meta; projection=\"3d\", origin=0.0, marker=\"+\", kwargs...)\n\nPlot mesh cell centers from axis view projection. projection should be either \"3d\", \"x\", \"y\" or \"z\". origin is center of projection plane in the normal direction.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.readmeta-Tuple{AbstractString}","page":"Internal","title":"Vlasiator.readmeta","text":"readmeta(filename; verbose=false) -> MetaData\n\nReturn MetaData from a vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readparameter-Tuple{MetaData, Any}","page":"Internal","title":"Vlasiator.readparameter","text":"readparameter(meta, param)\n\nReturn the parameter value from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvariable","page":"Internal","title":"Vlasiator.readvariable","text":"readvariable(meta, var, sorted=true) -> Array\n\nReturn variable value from the vlsv file. By default sorted=true, which means that for DCCRG grid the variables are sorted by cell ID.\n\n\n\n\n\n","category":"function"},{"location":"internal/#Vlasiator.readvariable-Tuple{MetaData, AbstractString, Any}","page":"Internal","title":"Vlasiator.readvariable","text":"readvariable(meta, var, idIn) -> Array\n\nRead a variable in a collection of cells.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvariablemeta-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.readvariablemeta","text":"readvariablemeta(meta, var) -> VarInfo\n\nReturn VarInfo about var in the vlsv file linked to meta.           \n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvcells-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.readvcells","text":"readvcells(meta, cellid; pop=\"proton\")\n\nRead velocity cells from a spatial cell of ID cellid, and return a map of velocity cell ids and corresponding value.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.refinedata-NTuple{5, Any}","page":"Internal","title":"Vlasiator.refinedata","text":"refinedata(meta, idlist, data, maxreflevel, normal) -> Array\n\nGenerate scalar data on the finest refinement level given cellids idlist and variable data on the slice perpendicular to normal.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.showdimension-Tuple{MetaData}","page":"Internal","title":"Vlasiator.showdimension","text":"showdimension(meta) -> Int\n\nReturn the dimension of data.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.showvariables-Tuple{MetaData}","page":"Internal","title":"Vlasiator.showvariables","text":"showvariables(meta) -> Vector{String}\n\nReturn all variable names in the VLSV file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#PyPlot-helpers","page":"Internal","title":"PyPlot helpers","text":"","category":"section"},{"location":"internal/","page":"Internal","title":"Internal","text":"Modules = [Vlasiator]\nPages   = [\"plot/pyplot.jl\"]","category":"page"},{"location":"internal/#Vlasiator.PlotArgs","page":"Internal","title":"Vlasiator.PlotArgs","text":"Plotting arguments.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.plot_prep2d-Tuple{Any, Any, Any, Any, Vlasiator.AxisUnit}","page":"Internal","title":"Vlasiator.plot_prep2d","text":"Generate axis and data for 2D plotting.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.set_args-Tuple{Any, Any, Vlasiator.AxisUnit, Vlasiator.ColorScale}","page":"Internal","title":"Vlasiator.set_args","text":"Set plot-related arguments.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.set_colorbar-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.set_colorbar","text":"Set colorbar norm and ticks.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.set_plot-NTuple{5, Any}","page":"Internal","title":"Vlasiator.set_plot","text":"Configure customized plot.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Private-APIs","page":"Internal","title":"Private APIs","text":"","category":"section"},{"location":"internal/","page":"Internal","title":"Internal","text":"Modules = [Vlasiator]\nPublic = false","category":"page"},{"location":"internal/#Vlasiator.AxisUnit","page":"Internal","title":"Vlasiator.AxisUnit","text":"Axis unit type. Currently supported: SI, RE.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.ColorScale","page":"Internal","title":"Vlasiator.ColorScale","text":"Color scales type for 2D plots. Currently supported: Log, Linear.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.MeshInfo","page":"Internal","title":"Vlasiator.MeshInfo","text":"Mesh size information.\n\n\n\n\n\n","category":"type"},{"location":"internal/#Vlasiator.getObjInfo-NTuple{5, Any}","page":"Internal","title":"Vlasiator.getObjInfo","text":"Return size and type information for the object.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getRotationB-Tuple{Any}","page":"Internal","title":"Vlasiator.getRotationB","text":"getRotationB(B) -> SMatrix\n\nObtain a rotation matrix with each column being a unit vector which is parallel (v3) and perpendicular (v1,v2) to the magnetic field B. The two perpendicular directions are chosen based on the reference vector of z-axis in the Cartesian coordinates.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.get_rotation_matrix-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.get_rotation_matrix","text":"get_rotation_matrix(vector, angle)\n\nCreates a rotation matrix that rotates around a unit vector by an angle in radians. References: https://en.wikipedia.org/wiki/Rodrigues'rotationformula https://en.wikipedia.org/wiki/Rotationmatrix#Rotationmatrixfromaxisandangle\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getfooter-Tuple{Any}","page":"Internal","title":"Vlasiator.getfooter","text":"Return the xml footer of vlsv.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.getindexes-NTuple{5, Any}","page":"Internal","title":"Vlasiator.getindexes","text":"Compute every cell id's x, y and z indexes on the given refinement level.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readmesh-NTuple{4, Any}","page":"Internal","title":"Vlasiator.readmesh","text":"Return mesh related variable.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.readvector-NTuple{4, Any}","page":"Internal","title":"Vlasiator.readvector","text":"Return vector data from vlsv file.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateTensorToVectorZ!-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.rotateTensorToVectorZ!","text":"rotateTensorToVector(tensor, vector)\n\nRotates tensor with a rotation matrix that aligns z-axis with vector.\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateWithB!-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.rotateWithB!","text":"rotateWithB!(T, B)\n\nRotate the tensor T with the 3rd direction aligned with B. See also: rotateWithB\n\n\n\n\n\n","category":"method"},{"location":"internal/#Vlasiator.rotateWithB-Tuple{Any, Any}","page":"Internal","title":"Vlasiator.rotateWithB","text":"rotateWithB(T, B) -> Matrix\n\nRotate the tensor T with the 3rd direction aligned with B. See also: rotateWithB!\n\n\n\n\n\n","category":"method"},{"location":"log/#Log","page":"Log","title":"Log","text":"","category":"section"},{"location":"log/#Performance","page":"Log","title":"Performance","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"The VLSV loader inherits the basic structure from Analysator and is redesigned for performance.","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Besides the language difference in speed, one of the key decisions in boosting performance is to avoid the usage of dictionary with integer keys as much as possible.","category":"page"},{"location":"log/#Benchmarks","page":"Log","title":"Benchmarks","text":"","category":"section"},{"location":"log/","page":"Log","title":"Log","text":"Initial tests on reading variables from sample VLSV files: ","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"DCCRG grid","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"2MB tmean [μs]\nJulia 200\nPython 1000","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"50MB tmean [μs]\nJulia 400\nPython 1000","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Field solver grid[1]","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 13\nPython 45","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"[1]: The field solver grid is a regular Cartesian grid at the finest refinement level. Therefore the storage requirement for fsgrid variables are quite significant: with 16 GB memory it is barely enough to read fg_b once; it will go out of memory for the second time!","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"I don't know why using Analysator is slower (2.3GB file, 4.8s) than directly calling matplotlib functions (2.3GB file, 0.5s). Same thing for Julia costs 1.0s (first time ~8s including everything).","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Reading and plotting one 2d slice of proton density out of 3D AMR data:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"26GB tmean [s]\nJulia 0.35\nPython 1.7","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Virtual satellite tracking from 845 frames of 3D AMR data (26G per frame) on Vorna:","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"1 CPU tmean [m]\nJulia 11\nPython 125","category":"page"},{"location":"log/","page":"Log","title":"Log","text":"Note that the above timings are for a single CPU. With only one command added for multithreading, the Julia timings can be improved by n where n is the number of threads. For example, with 8 threads, Julia takes ~80s to finish.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = Vlasiator","category":"page"},{"location":"#Vlasiator.jl","page":"Home","title":"Vlasiator.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Data processing and analyzing tool for the numerical model for collisionless ion-kinetic plasma physics Vlasiator. This lightweight package is built upon its sister in Python Analysator, and is carefully designed for performance and ease of use. It becomes even more powerful when using together with other fantastic packages: combined with external packages like FieldTracer.jl and TestParticle.jl, it is possible to do all kinds of in-depth analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Vlasiator.jl contains the following features:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Reading VLSV format data.\nExtracting quantities from the simulation at a given point/line/cut.\nPlotting 2D cuts.\nPlotting phase space distributions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThis package is still young, so be careful for any future breaking changes!","category":"page"},{"location":"","page":"Home","title":"Home","text":"warning: Warning\nThis package mostly aims at supporting Vlasiator 5.0+. Older versions of Vlasiator has different naming standard for outputs, and is not guaranteed to work.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install it,","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can then get started with","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Vlasiator","category":"page"},{"location":"","page":"Home","title":"Home","text":"More usages can be found in the examples.","category":"page"}]
}
