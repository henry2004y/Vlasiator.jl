"""
    viz(meta, var; args)

Visualize Vlasiator output `var` in `meta` with various options:

## Keyword arguments

* `axisunit`   - unit of axis of type `AxisUnit`
* `colorscale` - scale of colormap of type `ColorScale`
* `normal`     - slice normal direction
* `vmin`       - minimum color value
* `vmax`       - maximum color value
* `comp`       - selection of vector components

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function viz end

"""
    viz!(meta, var; args)

Visualize Meshes.jl `object` in an existing scene with `options` forwarded to [`viz`](@ref).
"""
function viz! end

function vlheatmap end

function vlslice end

"""
    vlslices(meta::MetaVLSV, var; axisunit=SI, comp=0, origin=[0.0, 0.0, 0.0])

Three orthogonal slices of `var` from `meta`.
"""
function vlslices end

"""
    vdfvolume(meta, location; species="proton", unit=SI, flimit=-1.0, verbose=false)

Meshscatter plot of VDFs in 3D.
"""
function vdfvolume end

function vdfslice end

function vdfslices end