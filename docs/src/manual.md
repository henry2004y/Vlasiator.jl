# Manual

Here we demonstrate some basic usages of processing Vlasiator output. For complete description of the arguments, please refer to the [API](internal.md) documents or type `?function_name` to display help message in the REPL.

## Common physical constants

A bunch of physical constants are predefined in [Vlasiator.jl](https://github.com/henry2004y/Vlasiator.jl/blob/master/src/vlsv/vlsvvariables.jl). To use them, you need to import explicitly, e.g. `using Vlasiator: RE` or prepend the module name like `Vlasiator.RE`.

| Physical constant | Value | Meaning |
|:---:|:--------------:|:-------------|
| qₑ | -1.60217662e-19 | electron charge, [C]             |
| mₑ | 9.10938356e-31  | electron mass, [kg]              |
| qᵢ | 1.60217662e-19  | proton mass, [C]                 |
| mᵢ | 1.673557546e-27 | proton mass, [kg]                |
| c  | 299792458.      | speed of light, [m/s]            |
| μ₀ | 4π*1e-7         | Vacuum permeability, [H/m]       |
| ϵ₀ | 1/(c^2*μ₀)      | Vacuum permittivity, [F/m]       |
| kB | 1.38064852e-23  | Boltzmann constant, [m²kg/(s²K)] |
| RE | 6.371e6         | Earth radius, [m]                |

## Loading VLSV data

- Read meta data

```julia
file = "bulk.0000004.vlsv"
meta = load(file)
```

The VLSV meta data contains information of file name, variable names, ordinary cell ID list, mesh sizes, species, and velocity cell structures. It is often the first argument for methods defined in Vlasiator.jl.

- Read parameter

For convenience we support the do-block syntax that automatically closes the file stream.

```julia
t = load(file) do meta
   readparameter(meta, "time")
end
```

- Read variable meta data

```julia
readvariablemeta(meta, "proton/vg_rho")
```

A list of utility functions has been implemented for checking variable status. See [here](internal.md#Vlasiator.hasname-Tuple{Any, Any, Any}) for the full list.

- Read variable

```julia
data = meta["proton/vg_rho"]
# Or equivalently
data = readvariable(meta, "proton/vg_rho")
```

The variable reading is designed for cells, which takes cell ID(s) as the 3rd argument if specified. The same interface works for both [DCCRG](https://github.com/fmihpc/dccrg) grid (for storing cell centered quantities like plasma moments) and [FS](https://github.com/fmihpc/fsgrid) grid (for storing field solver related quantities on a uniform high resolution mesh) variables. By default the returned DCCRG variable array is sorted by cell IDs. If in any case you want the original unsorted version as being stored in the file, use `readvariable(meta, var, false)`.

- Get variable at a given location

```julia
loc = [2.0, 0.0, 0.0]
id = getcell(meta, loc)
readvariable(meta, "proton/vg_rho", id)
```

- Get variable along a line between two points

```julia
using Vlasiator: RE # Earth radii
point1 = [12RE, 0, 0]
point2 = [15RE, 0, 0]
cellids, distances, coords = getcellinline(meta, point1, point2)
var_extract = readvariable(meta, "VA", cellids)
```

- Extract variable at a static cell ID from a sequence of files under the same grid

```julia
extractsat(files, var, id)
```

- Downsample field solver variable to DCCRG grid

```julia
Vlasiator.downsample_fg(meta, "fg_e")
```

- Upsample DCCRG variable to field solver grid

```julia
data = Vlasiator.read_variable_as_fg(meta, "proton/vg_rho")
# Equivalent to the above, but faster
data = let
   tmp = Vlasiator.fillmesh(meta, ["proton/vg_rho"]; maxamronly=true)[1][1][1]
   reshape(tmp, size(tmp)[2:end])
end
```

This is useful when corresponding DCCRG variables are not saved, or a uniform mesh is required for further analysis.

- Compare VLSV files

One may want to check if two vlsv files are identical. This is tricky because

1. the structure of VLSV format does not guarantee the writing order in parallel processing;
2. numerical error accumulates with floating point representation, especially with fastmath option.

The `issame` method does not check quantities that are related to the MPI writing sequence: for some reasons, even file sizes may vary depending on the number of MPI processes!

```julia
issame(file1, file2)
```

There is an optional third argument to `issame` for setting the relative difference tolerance, with default being 1e-4.
In practice relative difference works better for "large" numbers, and absolute difference works better for "small" numbers.

## Computing derived quantities

Vlasiator.jl is capable of computing plasma moments and some predefined derived quantities and saving them directly into VLSV files.
To avoid confusion about variable names, the conventions are

- raw quantities are all lowercases;
- all predefined derived variable names start with a capital letter;
- exceptions are for aliases (e.g. `"n"`).

To obtain a derived quantity, use either a key of string or symbol,

```julia
beta = meta["Beta"]
VA = meta[:VA]
```

Here is a full list of available quantities[^1]:

| Derived variable name | Meaning                          | Required variable[^2] |
|-----------------------|----------------------------------|-----------------------|
| Bmag                  | magnetic field magnitude         | vg\_b\_vol            |
| Emag                  | electric field magnitude         | vg\_e\_vol            |
| Vmag                  | bulk speed                       | vg\_v                 |
| Bhat                  | unit magnetic field              | vg\_b\_vol            |
| VS                    | sound speed                      | vg\_ptensor\_diagonal; vg\_rho |
| VA                    | Alfvén speed                     | vg\_rho; Bmag         |
| MA                    | Alfvén Mach number               | Vmag; VA              |
| MS                    | Sonic Mach number                | Vmag; VS              |
| Epar                  | $\mathbf{E}_\parallel$           | vg\_e\_vol; Bhat      |
| Eperp                 | $\mathbf{E}_\perp$               | vg\_e\_vol; Bhat      |
| Vpar                  | bulk velocity $\parallel\mathbf{B}$| vg\_v; vg\_b\_vol   |
| Vperp                 | bulk velocity $\perp \mathbf{B}$ | vg\_v; vg\_b\_vol     |
| Vth                   | proton thermal velocity          | P; vg\_rho            |
| P                     | scalar thermal pressure          | vg\_ptensor\_diagonal |
| Ppar                  | pressure $\parallel\mathbf{B}$   | vg\_ptensor\_diagonal; vg\_b\_vol |
| Pperp                 | pressure $\perp \mathbf{B}$      | vg\_ptensor\_offdiagonal; vg\_b\_vol |
| T                     | scalar temperature               | P; vg\_rho            |
| Tpar                  | temperature $\parallel\mathbf{B}$| vg\_rho; vg\_ptensor\_diagonal; vg\_b\_vol |
| Tperp                 | temperature $\perp \mathbf{B}$   | vg\_rho; vg\_ptensor\_offdiagonal; vg\_b\_vol |
| Tanisotropy           | $T_\perp / T_\parallel$          | Tpar; Tperp           |
| J                     | current density                  | vg\_b\_vol            |
| Jpar                  | $j_\parallel$                    | vg\_b\_vol            |
| Jperp                 | $j_\perp$                        | vg\_b\_vol            |
| Protated              | pressure tensor with $\widehat{z} \parallel \mathbf{B}$ | vg\_b\_vol; vg\_ptensor\_diagonal; vg\_ptensor\_offdiagonal |
| Panisotropy           | $P_\perp / P_\parallel$          | ptensor; B            |
| Pram                  | dynamic ram pressure             | vg\_rho; Vmag         |
| Pb                    | magnetic pressure                | vg\_b\_vol            |
| Poynting              | Poynting flux                    | E; B                  |
| Beta                  | plasma $\beta$, $P / P_B$        | P; vg\_b\_vol         |
| BetaStar              | modified $\beta$, $(P+P_{ram})/P_B$| P; Pram; vg\_b\_vol |
| IonInertial           | proton inertial length           | vg\_rho               | 
| Larmor                | proton Larmor radius             | Vth; Bmag             |
| Gyroperiod            | proton gyroperiod                | Bmag                  |
| PlasmaPeriod          | plasma oscillation period        | vg\_rho               |
| Gyrofrequency         | proton gyro-frequency            | Bmag                  |
| Omegap                | plasma frequency (proton)        | vg\_rho               |
| MagneticTension       | [magnetic tension force](https://en.wikipedia.org/wiki/Magnetic_tension) | vg\_b\_vol |
| n                     | proton number density            | vg\_rho               |

which can also be found as keys of dictionary in [vlsvvariables.jl](https://github.com/henry2004y/Vlasiator.jl/blob/master/src/vlsv/vlsvvariables.jl).

[^1]: For species specific variables, you need to add the species name at the front, separated by a slash. For example, the proton bulk velocity is a string `proton/vg_v`.
[^2]: If a required variable exists in the VLSV file, we try to use it directly instead of calculating from other variables. The interpolated FS grid variables onto DCCRG grid are preferred over original FS grid variables.

!!! note
    In Vlasiator, the cells inside the inner boundary (which is usually a sphere/circle) are filled with zero density values. This is then used to identify the inner boundary for all other quantities. Therefore, if you are manipulating directly on data, make sure that the nonsense values inside the inner boundary are excluded. One way to do this can be found in [vlsvvariables.jl](https://github.com/henry2004y/Vlasiator.jl/blob/master/src/vlsv/vlsvvariables.jl).

### Velocity space moments

We can also calculate plasma moments from the saved VLSV velocity space distributions.

```julia
# VDF cell indexes and values, with sparsity
vcellids, vcellf = readvcells(meta, cellid; species="proton")

getdensity(meta, vcellf)

getvelocity(meta, vcellids, vcellf)
# pressure tensor components Pxx, Pyy, Pzz, Pyz, Pzx, Pxy
getpressure(meta, vcellids, vcellf)
# heat flux components qⱼⱼᵢ
getheatfluxvector(meta, vcellids, vcellf)
```

To obtain the original ordering of velocity cells,

```julia
vcellids_original = Vlasiator.reorder(meta.meshes["proton"], vcellids)
```

Some useful quantities like non-Maxwellianity may be of interest. Currently we have implemented a monitor quantity named "Maxwellianity", which is defined as ``-ln \big[ 1/(2n) \int |f(v) - g(v)| dv \big]``, where n is the density, f(vᵢ) is the actual VDF value at velocity cell i, and g(vᵢ) is the analytical Maxwellian (or strictly speaking, normal) distribution with the same density, bulk velocity and scalar pressure as f.

```julia
getmaxwellianity(meta, vcellids, vcellf)
```

The value ranges from [0, +∞], with 0 meaning not Maxwellian-distributed at all, and +∞ a perfect Maxwellian distribution.

Sometimes it may be useful to recover the full 3D array of VDFs:

```julia
f = Vlasiator.reconstruct(meta.meshes["proton"], vcellids, vcellf)
```

However, usually in practice there would be only about 1% nonzero values. The moments and maxwellianity calculations above all have an alternative form of using reconstructed VDFs as inputs.

## Plotting

Vlasiator.jl does not include any plotting library as explicit dependency, but it offers plotting recipes/wrappers once the target plotting package is used.

Currently PyPlot provides the most complete and fine-tuned plotting capabilities.
Plots is a collection of plotting libraries with a uniform frontend, but it lacks detailed supports and consistent APIs.
Makie, a native Julia plotting library, is also supported via [VlasiatorMakie.jl](https://github.com/henry2004y/VlasiatorMakie.jl). Without generating an system image from [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl), it would take ~40s for the first plot on Julia 1.8. However, Makie has made nice progress in layouts, widgets, docs, and all the tiny things, which makes it a strong candidate for the de facto plotting library in the future.

More examples of customized plots can be found in the [repo](https://github.com/henry2004y/Vlasiator.jl/tree/master/src/examples).

### PyPlot Backend

To trigger Matplotlib plotting, `using PyPlot`.
All the functions with identical names as in Matplotlib accept all possible keyword arguments supported by their Matplotlib counterparts, e.g. font width, font size, colormap, etc. Users are encouraged to read the [Matplotlib](https://matplotlib.org/) documentation for details.

!!! warning
    The method call to certain axes is not dispatched, e.g. `ax.plot`; as an alternative, one needs to pass `ax` as the third argument to the functions, e.g. `plot(meta, "rho", ax)`. See [Matplotlib's two interfaces](https://aaltoscicomp.github.io/python-for-scicomp/data-visualization/#matplotlib-has-two-different-interfaces) for the history of the interfaces.

- Scalar colored contour from 2D simulation

```julia
pcolormesh(meta, "rho")
```

- Vector z-component colored contour from 2D simulation in a manually set range

```julia
pcolormesh(meta, "rho", comp=:z, colorscale=Log, axisunit=EARTH, vmin=1e6, vmax=2e6)
```

- Vz colored contour from 2D simulation with prescribed colormap

```julia
pcolormesh(meta, "proton/vg_v", comp=:z, colorscale=Linear, cmap=matplotlib.cm.RdBu_r)
```

- Derived quantity colored contour from 2D simulation (as long as the input variable is in the predefined dictionary)

```julia
pcolormesh(meta, "b", comp=:z, colorscale=Linear, axisunit=SI)
```

- Streamline from 2D simulation

```julia
streamplot(meta, "rho_v", comp="xy")
```

- Quiver from 2D simulation

```julia
quiver(meta, "rho_v", comp="xy")
```

The `comp` option is used to specify the two vector components.

You can choose to use a linear/log/symlog color scale by setting keyword `colorscale` to `Linear`, `Log`, or `SymLog`, plot vector components by setting keyword `op` to `:x`, `:y`, `:z`, `1`, `2`, `3`, `0` or `:mag`, and set `axisunit` to `EARTH` or `SI` etc.

- Mesh denoted by cell centers

```julia
plotmesh(meta; projection="z", color="w")
```

- Cut slice colored contour from 3D simulation

```julia
pcolormesh(meta, "proton/vg_rho", normal=:y, origin=0.0)
```

- Velocity distribution slice plot near a given spatial location `coordinates = [0.0, 0.0, 0.0]`

```julia
vdfslice(meta, coordinates)
```

- Extracted quantity line plot

```julia
rho_extract = vec(rho_extract)
loc = range(x1, x2, length=length(rho_extract))
plot(loc, rho_extract)
```

- Quick interactive REPL-based function for data inspection

```julia
pui(meta)
```

Or pass filename directly like `pui(file)`.

!!! note
    This is an experimental feature. We plan to have GUI-based plotting support in the future.

For a full list available optional arguments, please refer to the [doc for each method](internal.md#Public-APIs)

### Plots Backend

To trigger Plots.jl plotting, `using Plots`. This backend supports all available attributes provided by [Plots.jl](http://docs.juliaplots.org/latest/). By default it uses [GR](https://gr-framework.org/), but several other plotting libraries are also supported.

- Scaler colored contour from 2D simulation

```julia
heatmap(meta, var, aspect_ratio=:equal, c=:turbo)
```

- Scaler colored contour with lines from 2D simulation

```julia
contourf(meta, var)
```

- VDF projected slice in a normal direction

```julia
vdfslice(meta, location)
```

The keyword arguments are the same as in the PyPlot shown in the [API](internal.md#Vlasiator.vdfslice).

### Makie Backend

A standalone package [VlasiatorMakie.jl](https://github.com/henry2004y/VlasiatorMakie.jl) is designed for plotting with Makie. To trigger Makie plotting with OpenGL, `using VlasiatorMakie, GLMakie`.
You can either use intrinsic Makie plotting methods like

```julia
lines(meta, var)   # 1D
heatmap(meta, var) # 2D
```

or use full recipes provided by VlasiatorMakie

```julia
vlheatmap(meta, var)
```

For quick inspection of data, we have

- 2D slices of 3D AMR data

```julia
vlslice(meta, var; normal=:x)
```

- Orthognal slices of 3D AMR data

```julia
vlslices(meta, var)
```

- 2D slice of VDFs at a spatial cell

```julia
vdfslice(meta, location)
```

- Orthognal slices of VDFs at a spatial cell

```julia
vdfslices(meta, location)
```

- 3D scatter of VDFs at a spatial cell

```julia
vdfvolume(meta, location)
```

The interactive plots are available through the OpenGL backend of Makie `GLMakie`. For noninteractive high fidelity plots, we can also use the Cairo backend of Makie `CairoMakie`. Other options can be found at [Makie Ecosystem](https://makie.juliaplots.org/stable/#makie_ecosystem).

## Appending to VLSV

We are able to compute derived quantities from an original VLSV file and generate a new VLSV output with new quantities included.

```julia
vmag = readvariable(meta, "Vmag", meta.cellid)
pa = readvariable(meta, "Panisotropy", meta.cellid)
vars = Vector{Tuple{VecOrMat, String, VarInfo}}(undef, 0)
# require LaTeXStrings.jl
push!(vars, (vmag, "vmag", VarInfo("m/s", L"$\mathrm{m}/mathrm{s}$", L"$V$", "")))
push!(vars, (pa, "panisotropy", VarInfo("", "", "", "")))

write_vlsv("bulk.vlsv", "bulk_new.vlsv", vars)
```

!!! note
    Writing new FsGrid variables is not supported. All quantities from the original file is maintained.

## Converting to VTK

We can convert VLSV files into VTK format. Since DCCRG is Cartesian based with uniform spacing, each level of mesh refinement corresponds to a VTK image file, and the cell refinement relationships are stored in `vtkGhostType` as well as the `vthb` file.

To convert a VLSV file into VTK,

```julia
write_vtk(file)
```

This function accepts either string of file names or `MetaVLSV`.

To see the full list of options, please refer to the documentation in [API Reference](internal.md). Example usage can be found [here](https://github.com/henry2004y/Vlasiator.jl/blob/master/examples/demo_convert2vti.jl).

!!! warning
    As of ParaView 5.9.1, there are [display issues](https://discourse.paraview.org/t/vthb-file-structure/7224) with `VTKOverlappingAMR`. However, we can read the generated image files directly. There is also an keyword argument for `write_vtk` called `maxamronly`: when it is set to `true`, then only the image file at the highest refinement level is generated.
    This part is experimental and subject to change in the future.

## Tracking log files

The runtime performance per iteration can be monitored through log files:

```julia
file = "logfile.txt"
timestamps, speed = readlog(file)
```

See a live example at [demo_log.jl](https://github.com/henry2004y/Vlasiator.jl/tree/master/examples/demo_log.jl).

## Examples

More [examples](https://github.com/henry2004y/Vlasiator.jl/tree/master/examples) are provided about

- Plotting with PyPlot
- Plotting with Plots
- Extracing variable along a line
- Field line tracing
- Simulation log file tracking
- Converting VLSV to VTK format
- Parallel post-processing
- Finding X-points and O-points in 2D reconnections

Feel free to check those out and try on your data!
