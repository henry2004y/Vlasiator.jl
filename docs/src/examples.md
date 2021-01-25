# Examples

## Loading VLSV data

- Read meta data
```
filename = "bulk.0000004.vlsv"
meta = read_meta(filename)
```

- Read variable
```
data = read_variable(meta, "proton/vg_rho")
```

- Get variable at a given location (This can be simplified even further later!)
```
loc = [2.0, 0.0, 0.0]
id = get_cellid(meta, loc)
read_variable_select(meta, "proton/vg_rho", id)
```

- Get variable along a line between two points
```
point1 = [12Re, 0, 0]
point2 = [15Re, 0, 0]
cellids, distances, coords = get_cell_in_line(meta, point1, point2)
```

## Plotting

- Scaler colored contour for 2D simulation
```
plot_pcolormesh(meta, "rho")
```

- Vector z component colored contour for 2D simulation
```
plot_pcolormesh(meta, "rho", op="z", islinear=false, axisunit="Re")
```

- Derived quantity colored contour for 2D simulation (as long as the input variable is in the predefined dictionary)
```
plot_pcolormesh(meta, "b", op="z", islinear=false, axisunit="Re")
```

- Streamline for 2D simulation
```
streamline(meta, "rho_v", comp="xy")
```

The `comp` option is used to specify the two vector components.

!!! note
    Currently there is limited support for derived variables. This will be expanded and changed later for ease of use!

You can choose to use linear/log color scale, plot vector components via e.g. `op="x"` or magnitude by default, and set unit via `axisunit="Re"` etc..

- Cut slice colored contour for 3D simulation
```
plot_colormap3dslice(meta, "proton/vg_rho", normal="y")
```

More examples of customized plots can be found in the [repo](https://github.com/henry2004y/Vlasiator.jl/tree/master/src/examples).