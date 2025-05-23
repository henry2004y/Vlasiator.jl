"Plotting functionalities from Matplotlib."
module VlasiatorPyPlot

using Reexport
@reexport using Vlasiator
using Vlasiator: AxisUnit, ColorScale, PlotArgs,
   set_args, prep2d, prep2dslice, get_axis, _fillinnerBC!, _getdim2d, set_lim, prep_vdf,
   SVector
@reexport using PyPlot

import PyPlot.PyCall: PyObject
using REPL.TerminalMenus # Command line UI

import Vlasiator: vdfslice
export plot, pcolormesh, streamplot, quiver, plotmesh, set_colorbar
export pui

function __init__()
   if VersionNumber(matplotlib.__version__) ≥ v"3.3"
      matplotlib.rc("image", cmap="turbo") # set default colormap
   end

   if VersionNumber(matplotlib.__version__) < v"3.5"
      matplotlib.rc("pcolor", shading="nearest") # newer version default "auto"
   end

   matplotlib.rc("font", size=15)
   matplotlib.rc("xtick", labelsize=10)
   matplotlib.rc("ytick", labelsize=10)
   matplotlib.rc("xtick.minor", visible=true)
   matplotlib.rc("ytick.minor", visible=true)
end

"""
    plot(meta::MetaVLSV, var::String, ax=nothing; comp=0, kwargs...)

Plot 1D `var` from data linked to VLSV `meta`. If `ax===nothing`, plot on the current axes.
The keyword `comp` is used to specify the component index of a vector, with 0 being the
magnitude. Any valid keyword argument for `plt.plot` is accepted.

# Keywords
- `comp::Int`: the component of a vector, chosen from any valid integers.
"""
function PyPlot.plot(meta::MetaVLSV, var::String, ax::Union{PyObject,Nothing}=nothing;
   comp::Int=0, kwargs...)
   ndims(meta) == 1 || error("plot only accepts 1D data!")

   datafull = readvariable(meta, var)

   if ndims(datafull) == 1
      data = datafull
   else
      if comp == 0
         data = sqrt.(sum(x->x^2, datafull, dims=1)) |> vec
      else
         data = datafull[comp,:]
      end
   end

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   isnothing(ax) && (ax = plt.gca())

   ax.plot(x, data; kwargs...)
end

"""
    streamplot(meta::MetaVLSV, var::String, ax=nothing; comp::String="xy", axisunit=EARTH,
       kwargs...)

Wrapper over Matplotlib's streamplot function.

# Keywords
- `comp::String`: a subset of "xyz" in any order.
- `axisunit::AxisUnit`: chosen from `EARTH` and `SI`.

Any valid keyword argument for `plt.streamplot` is accepted.
"""
function PyPlot.streamplot(meta::MetaVLSV, var::String, ax::Union{PyObject,Nothing}=nothing;
   comp::String="xy", axisunit::AxisUnit=EARTH, origin::AbstractFloat=0.0, kwargs...)

   x, y, v1, v2 = set_vector(meta, var, comp, axisunit, origin)

   isnothing(ax) && (ax = plt.gca())

   ax.streamplot(x, y, v1, v2; kwargs...)
end

"""
    quiver(meta::MetaVLSV, var::String, ax=nothing; comp="xy", axisunit=EARTH, stride=10,
       kwargs...)

Wrapper over Matplotlib's quiver function. If `ax===nothing`, plot on the current axes.

# Keywords
- `comp::String`: a subset of "xyz" in any order.
- `axisunit::AxisUnit`: chosen from `EARTH` and `SI`.
- `stride::Int`: arrow strides in number of cells.

Any valid keyword argument for `plt.quiver` is accepted.
"""
function PyPlot.quiver(meta::MetaVLSV, var::String, ax::Union{PyObject,Nothing}=nothing;
   comp::String="xy", axisunit::AxisUnit=EARTH, origin::AbstractFloat=0.0, stride::Int=10,
   kwargs...)
   x, y, v1, v2 = set_vector(meta, var, comp, axisunit, origin)

   isnothing(ax) && (ax = plt.gca())

   Xq,  Yq  = x[1:stride:end], y[1:stride:end]
   v1q, v2q = v1[1:stride:end, 1:stride:end], v2[1:stride:end, 1:stride:end]

   ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

function set_vector(meta::MetaVLSV, var::String, comp::String, axisunit::AxisUnit,
   origin::AbstractFloat = 0.0)
   @assert ndims(meta) == 2 "Vector plot only support 2D outputs!"
   (; ncells, maxamr, coordmin, coordmax) = meta
   if occursin("x", comp)
      v1_ = 1
      if occursin("y", comp)
         dir = 3
         v2_ = 2
         sizes = _getdim2d(ncells, maxamr, :z)
         plotrange = (coordmin[1], coordmax[1], coordmin[2], coordmax[2])
      else
         dir = 2
         v2_ = 3
         sizes = _getdim2d(ncells, maxamr, :y)
         plotrange = (coordmin[1], coordmax[1], coordmin[3], coordmax[3])
      end
   else
      dir = 1
      v1_, v2_ = 2, 3
      sizes = _getdim2d(ncells, maxamr, :x)
      plotrange = (coordmin[2], coordmax[2], coordmin[3], coordmax[3])
   end

   data = readvariable(meta, var)

   if !startswith(var, "fg_") # vlasov grid
      @assert ndims(data) == 2 && size(data,1) == 3 "Vector variable required!"
      if meta.maxamr == 0
         data = reshape(data, 3, sizes[1], sizes[2])
         v1 = data[v1_,:,:]
         v2 = data[v2_,:,:]
      else
         sliceoffset = origin - coordmin[dir]
         idlist, indexlist =
            getslicecell(meta, sliceoffset, dir, coordmin[dir], coordmax[dir])
         v2D = data[:,indexlist]
         v1 = @views refineslice(meta, idlist, v2D[v1_,:], dir)
         v2 = @views refineslice(meta, idlist, v2D[v2_,:], dir)
      end
   else # Field solver grid
      dataslice = dropdims(data, dims = tuple(findall(size(data) .== 1)...))
      v1 = dataslice[v1_,:,:]
      v2 = dataslice[v2_,:,:]
   end

   x, y = get_axis(axisunit, plotrange, sizes)

   x, y, v1', v2'
end

"""
    pcolormesh(meta::MetaVLSV, var::String, ax=nothing;
       comp=0, axisunit=EARTH, colorscale=Linear, vmin=-Inf, vmax=Inf,
       addcolorbar=true, extent=[-Inf, Inf, -Inf, Inf], kwargs...)

Plot a variable using pseudocolor from `meta` associated with 2D or 3D AMR grid.
If `ax` is provided, then it will plot on that axes.

# Keywords
- `comp::Union{Int, Symbol}`: the component of a vector, chosen from `:mag, :x, :y, :z, 0`
or valid integers for indexing.
- `axisunit::AxisUnit`: the unit of axis ∈ `EARTH, SI`.
- `colorscale::ColorScale`: `Linear`, `Log`, or `SymLog`.
- `vmin::Float64`: minimum data range. Set to maximum of data if not specified.
- `vmax::Float64`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.
- `extent::Vector{Float64}`: extent of axis ranges in the same unit as `axisunit`.
- `normal::Symbol`: normal direction for slice of 3D data, `:x`, `:y`, `:z`.
- `origin::Float64`: origin of plane slice of 3D data.

Any valid keyword argument for `plt.pcolormesh` is supported.

# Example

```julia
pcolormesh(meta, var)

pcolormesh(meta, var, axisunit=SI)

pcolormesh(data, var, colorscale=Log, extent=[0,1,0,2])
```
"""
function PyPlot.pcolormesh(meta::MetaVLSV, var::String,
   ax::Union{PyObject,Nothing}=nothing; comp::Union{Int, Symbol}=0,
   axisunit::AxisUnit=EARTH, colorscale::ColorScale=Linear, addcolorbar::Bool=true,
   vmin::Float64=-Inf, vmax::Float64=Inf, extent::Vector{Float64}=[-Inf, Inf, -Inf, Inf],
   kwargs...)

   isnothing(ax) && (ax = plt.gca())

   _plot2d(ax.pcolormesh, meta, var, ax;
      comp, axisunit, colorscale, addcolorbar, vmin, vmax, extent, kwargs...)
end

"""
    contour(meta::MetaVLSV, var::String, ax=nothing;
       comp=0, axisunit=EARTH, colorscale=Linear, vmin=-Inf, vmax=Inf,
       addcolorbar=true, extent=[-Inf, Inf, -Inf, Inf], kwargs...)

Plot contour lines of `var` from `meta` associated with 2D or 3D AMR grid.
If `ax` is provided, then it will plot on that axes. The arguments are the same as in
[`pcolormesh`](@ref).

Any valid keyword argument for `plt.contour` is supported.
"""
function PyPlot.contour(meta::MetaVLSV, var::String,
   ax::Union{PyObject,Nothing}=nothing; comp::Union{Int, Symbol}=0,
   axisunit::AxisUnit=EARTH, colorscale::ColorScale=Linear, addcolorbar::Bool=true,
   vmin::Float64=-Inf, vmax::Float64=Inf, extent::Vector{Float64}=[-Inf, Inf, -Inf, Inf],
   kwargs...)

   isnothing(ax) && (ax = plt.gca())

   _plot2d(ax.contour, meta, var, ax;
      comp, axisunit, colorscale, addcolorbar, vmin, vmax, extent, kwargs...)
end

"""
    contourf(meta::MetaVLSV, var::String, ax=nothing;
       comp=0, axisunit=EARTH, colorscale=Linear, vmin=-Inf, vmax=Inf,
       addcolorbar=true, extent=[-Inf, Inf, -Inf, Inf], kwargs...)

Plot colored contour lines of `var` from `meta` associated with 2D or 3D AMR grid.
If `ax` is provided, then it will plot on that axes. The arguments are the same as in
[`pcolormesh`](@ref). See also [`contour`](@ref).

Any valid keyword argument for `plt.contourf` is supported.
"""
function PyPlot.contourf(meta::MetaVLSV, var::String,
   ax::Union{PyObject,Nothing}=nothing; comp::Union{Int, Symbol}=0,
   axisunit::AxisUnit=EARTH, colorscale::ColorScale=Linear, addcolorbar::Bool=true,
   vmin::Float64=-Inf, vmax::Float64=Inf, extent::Vector{Float64}=[-Inf, Inf, -Inf, Inf],
   kwargs...)

   isnothing(ax) && (ax = plt.gca())

   _plot2d(ax.contourf, meta, var, ax;
      comp, axisunit, colorscale, addcolorbar, vmin, vmax, extent, kwargs...)
end

# General method for 2D plotting
function _plot2d(f::PyObject, meta::MetaVLSV, var::String,
   ax::Union{PyObject,Nothing}=nothing; comp::Union{Int, Symbol}=0,
   axisunit::AxisUnit=EARTH, colorscale::ColorScale=Linear, addcolorbar::Bool=true,
   vmin::Float64=-Inf, vmax::Float64=Inf, extent::Vector{Float64}=[-Inf, Inf, -Inf, Inf],
   kwargs...)

   if ndims(meta) == 3 || meta.maxamr > 0
      # check if origin and normal exist in kwargs
      normal = haskey(kwargs, :normal) ? values(kwargs).normal : :y
      origin = haskey(kwargs, :origin) ? values(kwargs).origin : 0.0
      kwargs = Base.structdiff(values(kwargs), (normal = normal, origin = origin))

      pArgs = set_args(meta, var, axisunit; normal, origin)
      data = prep2dslice(meta, var, normal, comp, pArgs)'
   else
      pArgs = set_args(meta, var, axisunit)
      data = prep2d(meta, var, comp)'
   end

   x1, x2 = get_axis(pArgs)

   if var in ("fg_b", "fg_e", "vg_b_vol", "vg_e_vol") || endswith(var, "vg_v")
      _fillinnerBC!(data, data)
   end

   norm, ticks = set_colorbar(colorscale, vmin, vmax, data)

   range1 = searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2])
   range2 = searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])

   if colorscale != SymLog
      if range1 == 1:pArgs.sizes[1] && range2 == 1:pArgs.sizes[2]
         c = f(x1, x2, data; norm, kwargs...)
      else
         c = f(x1[range1], x2[range2], data[range2, range1];
            norm, kwargs...)
      end
   else
      if range1 == 1:pArgs.sizes[1] && range2 == 1:pArgs.sizes[2]
         c = f(x1, x2, data; norm, cmap=matplotlib.cm.RdBu_r, kwargs...)
      else
         c = f(x1[range1], x2[range2], data[range2, range1];
            norm, cmap=matplotlib.cm.RdBu_r, kwargs...)
      end
   end

   set_plot(c, ax, pArgs, ticks, addcolorbar)

   return c
end

"""
    set_colorbar(colorscale::ColorScale, v1, v2, data=[1.0;;];
       linthresh=1.0, logstep=1) -> norm, ticks

Set colorbar norm and ticks in a given range `v1` to `v2` for `data` in `colorscale`.
Matplotlib's Colormap Normalization Section presents various kinds of normlizations beyond
linear, logarithmic and symmetric logarithmic, like centered, discrete, and two slope norm.
For fine-grained control, use the norm methods from `matplotlib.colors` directly.

# Example

```julia
norm = matplotlib.colors.CenteredNorm(); # after matplotlib v3.4
norm = matplotlib.colors.BoundaryNorm(boundaries=[0, 1], ncolors=2);
norm = matplotlib.colors.TwoSlopeNorm(vmin=-500., vcenter=0, vmax=4000);
norm = matplotlib.colors.SymLogNorm(;linthresh=1.0, linscale=0.01, vmin=-1e2, vmax=1e2, base=10)
```
There are also various options for the ticks available in `matplotlib.ticker`.
"""
function set_colorbar(colorscale::ColorScale, v1::AbstractFloat, v2::AbstractFloat,
   data::AbstractArray{<:Real, 2}=[1.0;;];
   linthresh::AbstractFloat=1.0, logstep::AbstractFloat=1.0, linscale::AbstractFloat=0.03)

   vmin, vmax = set_lim(v1, v2, data, colorscale)
   if colorscale == Linear
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      norm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
      ticks = matplotlib.ticker.LinearLocator(numticks=9)
   elseif colorscale == Log # logarithmic
      norm = matplotlib.colors.LogNorm(;vmin, vmax)
      ticks = matplotlib.ticker.LogLocator(base=10, subs=collect(0:9))
   else # symmetric log
      logthresh = floor(Int, log10(linthresh))
      minlog = ceil(Int, log10(-vmin))
      maxlog = ceil(Int, log10(vmax))

      norm = matplotlib.colors.SymLogNorm(;linthresh, linscale, vmin, vmax, base=10)
      ticks = [ [-(10.0^x) for x in minlog:-logstep:logthresh]..., 0.0,
         [10.0^x for x in logthresh:logstep:maxlog]..., ]
   end

   norm, ticks
end

"Configure customized plot."
function set_plot(c::PyObject, ax::PyObject, pArgs::PlotArgs, ticks, addcolorbar::Bool)
   (;str_title, strx, stry, cb_title) = pArgs

   if addcolorbar
      cb = plt.colorbar(c; ax, ticks, fraction=0.04, pad=0.02)
      !isempty(cb_title) && cb.ax.set_ylabel(cb_title)
      cb.ax.tick_params(direction="in")
   end

   ax.set_title(str_title, fontweight="bold")
   ax.set_xlabel(strx)
   ax.set_ylabel(stry)
   ax.set_aspect("equal")

   # Set border line widths
   for loc in ("left", "bottom", "right", "top")
      VersionNumber(matplotlib.__version__) < v"3.4" ?
         ax.spines[loc].set_linewidth(2.0) :
         getproperty(ax.spines, loc).set_linewidth(2.0)
   end

   ax.xaxis.set_tick_params(width=2.0, length=3)
   ax.yaxis.set_tick_params(width=2.0, length=3)
   return
end

"""
    vdfslice(meta::MetaVLSV, location::AbstractVector, ax=nothing; kwargs...)

Plot the 2D slice cut of phase space distribution function at `location` within velocity
range `limits`. If `ax===nothing`, plot on the current active axes.

# Keywords

- `unit::AxisUnit`: location unit in `SI`, `EARTH`.
- `unitv::String`: velocity unit in ("km/s", "m/s").
- `limits::Vector{Real}`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype::Symbol`: slice type from `:xy`, `:xz`, `:yz`, `:bperp`, `:bpar1`, `:bpar2`.
- `center::Symbol`: setting the reference frame from `:bulk`, `:peak`.
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected!
- `vmin, vmax::AbstractFloat`: minimum and maximum VDF values for plotting.
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
- `flimit::AbstractFloat`: minimum VDF threshold for plotting.

Any valid keyword argument for `hist2d` is supported.
"""
function vdfslice(meta::MetaVLSV, location::AbstractVector,
   ax::Union{PyObject,Nothing}=nothing; limits::Vector{Float64}=[-Inf, Inf, -Inf, Inf],
   verbose::Bool=false, species::String="proton", unit::AxisUnit=SI, unitv::String="km/s",
   vmin::AbstractFloat=-Inf, vmax::AbstractFloat=Inf, slicetype::Symbol=:default,
   vslicethick::AbstractFloat=0.0, center::Symbol=:nothing, weight::Symbol=:particle,
   flimit::AbstractFloat=-1.0, kwargs...)

   v1, v2, r1, r2, weights, strx, stry, str_title =
      prep_vdf(meta, location;
         species, unit, unitv, slicetype, vslicethick, center, weight, flimit, verbose)

   isinf(vmin) && (vmin = minimum(weights))
   isinf(vmax) && (vmax = maximum(weights))

   verbose && @info "Active f range is $vmin, $vmax"

   isnothing(ax) && (ax = plt.gca())

   norm = matplotlib.colors.LogNorm(vmin, vmax)

   h = ax.hist2d(v1, v2; bins=(r1, r2), weights, norm, shading="flat")

   ax.set_title(str_title, fontweight="bold")
   ax.set_xlabel(strx)
   ax.set_ylabel(stry)
   ax.set_aspect("equal")
   ax.grid(color="grey", linestyle="-")
   ax.tick_params(direction="in")

   cb = plt.colorbar(h[4]; ax, fraction=0.04, pad=0.02)
   cb.ax.tick_params(which="both", direction="in")
   cb_title = cb.ax.set_ylabel("f(v)")

   if slicetype in (:xy, :xz, :yz)
      # Draw vector of magnetic field direction
   end

   h[4] # h[1] is 2D data, h[2] is x axis, h[3] is y axis
end

"""
    plotmesh(meta::MetaVLSV; projection::String="3d", origin::AbstractFloat=0.0,
       marker::String="+", kwargs...)

Plot mesh cell centers from axis view `projection`. `projection` should be either "3d", "x",
"y" or "z". `origin` is center of projection plane in the normal direction.
"""
function plotmesh(meta::MetaVLSV, ax::Union{PyObject,Nothing}=nothing;
   projection::String="3d", origin::AbstractFloat=0.0, marker::String="+", kwargs...)
   (;coordmin, coordmax) = meta
   if projection in ("x", "y", "z")
      dirp, dir1, dir2 =
         if projection == "x"
            1, 2, 3
         elseif projection == "y"
            2, 1, 3
         else # "z"
            3, 1, 2
         end
      sliceoffset = origin - coordmin[dirp]
      ids, _ = getslicecell(meta, sliceoffset, dirp, coordmin[dirp], coordmax[dirp])
   else # 3D
      ids = Vlasiator.getcellid(meta.fid, meta.nodeVLSV.var)
   end

   centers = [zeros(SVector{3, Float32}) for _ in ids]
   for (i, id) in enumerate(ids)
      @inbounds centers[i] = getcellcoordinates(meta, id)
   end

   isnothing(ax) && (ax = plt.gca())

   if projection in ("x", "y", "z")
      x1 = getindex.(centers, dir1)
      x2 = getindex.(centers, dir2)
      s = ax.scatter(x1, x2; marker, kwargs...)
   else
      x1 = getindex.(centers, 1)
      x2 = getindex.(centers, 2)
      x3 = getindex.(centers, 3)

      if ax.name == "3d"
         s = ax.scatter(x1, x2, x3; marker, kwargs...)
      else
         @error "Keyword projection set to \"3d\". Please use 3d projection axis!"
      end
   end
   s
end

"""
    pui(meta::MetaVLSV; suppress_output=false)
    pui(file::AbstractString)

Quick plotting via command line interactive selections.
"""
function pui(meta::MetaVLSV; suppress_output::Bool=false)
   menu = RadioMenu(meta.variable; charset=:ascii)

   println("Choose variable to plot:")
   var_ = request(menu; suppress_output)

   var_ == -1 && error("Variable selection canceled.")

   comp = 0

   if meta.variable[var_] in ("fg_b", "fg_e", "vg_b_vol", "vg_e_vol") ||
      endswith(meta.variable[var_], "vg_v") ||
      occursin("vg_ptensor", meta.variable[var_])

      menu = RadioMenu(["x","y","z","mag"]; charset=:ascii)

      comp = request("Choose vector component to plot:", menu; suppress_output)

      comp == -1 && error("Vector component selection canceled.")

      comp == 4 && (comp = 0)
   end

   if ndims(meta) != 1
      pcolormesh(meta, meta.variable[var_]; comp)
   else # 1D
      plot(meta, meta.variable[var_]; comp)
   end

   return
end

pui(file::AbstractString) = file |> load |> plot_ui

end