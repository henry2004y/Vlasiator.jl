# Plotting functionalities from Matplotlib.

using PyPlot
using REPL.TerminalMenus # Command line UI

export plot, pcolormesh, pcolormeshslice, vdfslice, streamplot, quiver, plotmesh
export pui

@static if matplotlib.__version__ ≥ "3.3"
   matplotlib.rc("image", cmap="turbo") # set default colormap
end

@static if matplotlib.__version__ < "3.5"
   matplotlib.rc("pcolor", shading="nearest") # newer version default "auto"
end

matplotlib.rc("font", size=14)
matplotlib.rc("xtick", labelsize=10)
matplotlib.rc("ytick", labelsize=10)


"""
    plot(meta, var, ax=nothing; comp=0, kwargs)

Plot `var` from `meta` of 1D VLSV data. If `ax===nothing`, plot on the current active axes.
The keyword `comp` in (0, 1, 2, 3) is used to specify the component of a vector.
"""
function PyPlot.plot(meta::MetaVLSV, var, ax=nothing; comp=0, kwargs...)
   ndims(meta) == 1 || error("plot only accepts 1D data!")

   datafull = readvariable(meta, var)

   if comp == 0
      data = datafull
   else
      data = datafull[comp,:]
   end

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   isnothing(ax) && (ax = plt.gca())

   ax.plot(x, data; kwargs...)
end

"""
    streamplot(meta, var, ax=nothing; comp="xy", axisunit=EARTH, kwargs...)

Wrapper over Matplotlib's streamplot function.
The keyword arguments can be any valid Matplotlib arguments into streamplot.

# Optional arguments
- `comp`: a subset of "xyz" in any order.
- `axisunit`: chosen from `EARTH` and `SI`.
"""
function PyPlot.streamplot(meta::MetaVLSV, var::AbstractString, ax=nothing;
   comp="xy", axisunit=EARTH, kwargs...)

   X, Y, v1, v2 = set_vector(meta, var, comp, axisunit)

   isnothing(ax) && (ax = plt.gca())

   streamplot(X, Y, v1, v2; kwargs...)
end

"""
    quiver(meta, var, ax=nothing; comp="xy", axisunit=EARTH, stride=10, kwargs...)

Wrapper over Matplotlib's quiver function. If `ax===nothing`, plot on the current active
axes.
The keyword arguments can be any valid Matplotlib arguments into quiver.

# Optional arguments
- `comp`: a subset of "xyz" in any order.
- `axisunit`: chosen from `EARTH` and `SI`.
- `stride::Integer`: arrow strides in number of cells.
"""
function PyPlot.quiver(meta::MetaVLSV, var::AbstractString, ax=nothing;
   comp="xy", axisunit::AxisUnit=EARTH, stride::Integer=10, kwargs...)

   X, Y, v1, v2 = set_vector(meta, var, comp, axisunit)

   isnothing(ax) && (ax = plt.gca())

   Xq,  Yq  = X[1:stride:end, 1:stride:end],  Y[1:stride:end, 1:stride:end]
   v1q, v2q = v1[1:stride:end, 1:stride:end], v2[1:stride:end, 1:stride:end]

   ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

function set_vector(meta::MetaVLSV, var, comp, axisunit::AxisUnit)
   (;ncells, coordmin, coordmax) = meta
   if occursin("x", comp)
      v1_ = 1
      if occursin("y", comp)
         v2_ = 2
         sizes = (ncells[1], ncells[2])
         plotrange = (coordmin[1], coordmax[1], coordmin[2], coordmax[2])
      else
         v2_ = 3
         sizes = (ncells[1], ncells[3])
         plotrange = (coordmin[1], coordmax[1], coordmin[3], coordmax[3])
      end
   else
      v1_, v2_ = 2, 3
      sizes = (ncells[2], ncells[3])
      plotrange = (coordmin[2], coordmax[2], coordmin[3], coordmax[3])
   end

   data = readvariable(meta, var)

   if !startswith(var, "fg_") # vlasov grid
      @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required!"
      data = reshape(data, 3, sizes[1], sizes[2])
   end

   v1 = data[v1_,:,:]
   v2 = data[v2_,:,:]

   x, y = get_axis(axisunit, plotrange, sizes)

   # meshgrid: note the array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   X, Y, v1', v2'
end

"""
    pcolormesh(meta::MetaVLSV, var::AbstractString, ax=nothing;
       comp=0, axisunit=EARTH, colorscale=Linear, vmin=-Inf32, vmax=Inf32,
       addcolorbar=true, extent=[-Inf32, Inf32, -Inf32, Inf32], kwargs...)

Plot a variable using pseudocolor from 2D VLSV data.
If `ax` is provided, then it will plot on that axes.
If 3D or AMR grid detected, it will pass arguments to [`pcolormeshslice`](@ref).

# Optional arguments
- `comp::Tuple{Int64, Symbol}`: the component of a vector, chosen from
`:mag, :x, :y, :z, 0, 1, 2, 3`.
- `axisunit::AxisUnit`: the unit of axis ∈ `EARTH, SI`.
- `colorscale::ColorScale`: `Linear`, `Log`, or `SymLog`.
- `vmin::Float`: minimum data range. Set to maximum of data if not specified.
- `vmax::Float`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.
- `extent::Vector`: extent of axis ranges for plotting in the same unit as `axisunit`.
- `normal::Symbol`: normal direction for slice of 3D data, `:x`, `:y`, `:z`.
- `origin::Float`: origin of plane slice of 3D data.

`pcolormesh(meta, var)`

`pcolormesh(meta, var, axisunit=SI)`

`pcolormesh(data, var, colorscale=Log, extent=[0,1,0,2])`
"""
function PyPlot.pcolormesh(meta::MetaVLSV, var::AbstractString, ax=nothing;
   comp=0, axisunit::AxisUnit=EARTH, colorscale::ColorScale=Linear, addcolorbar=true,
   vmin=-Inf32, vmax=Inf32, extent=[-Inf32, Inf32, -Inf32, Inf32], kwargs...)

   if ndims(meta) == 3 || meta.maxamr > 0
      # check if origin and normal exist in kwargs
      normal = haskey(kwargs, :normal) ? kwargs.data.normal : :y
      origin = haskey(kwargs, :origin) ? kwargs.data.origin : 0.0
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

   range1, range2 =
      if axisunit == EARTH
         searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
         searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])
      else
         searchsortedfirst(x1, extent[1]):searchsortedlast(x1, extent[2]),
         searchsortedfirst(x2, extent[3]):searchsortedlast(x2, extent[4])
      end

   isnothing(ax) && (ax = plt.gca())

   if colorscale != SymLog
      c = ax.pcolormesh(x1[range1], x2[range2], data[range2, range1];
         norm, kwargs...)
   else
      c = ax.pcolormesh(x1[range1], x2[range2], data[range2, range1];
         norm, cmap=matplotlib.cm.RdBu_r, kwargs...)
   end

   set_plot(c, ax, pArgs, ticks, addcolorbar)

   return c
end

"""
    set_colorbar(colorscale::ColorScale, v1, v2, data=[1.0]) -> norm, ticks

Set colorbar norm and ticks in a given range `v1` to `v2` for `data` in `colorscale`.
Matplotlib's Colormap Normalization Section presents various kinds of normlizations beyond
linear, logarithmic and symmetric logarithmic, like centered, discrete, and two slope norm.
For fine-grain control, it is suggested to use the norm methods from `matplotlib.colors`.
For instance,
```
julia> norm = matplotlib.colors.CenteredNorm(); # after matplotlib v3.4
julia> norm = matplotlib.colors.BoundaryNorm(boundaries=[0, 1], ncolors=2);
julia> norm = matplotlib.colors.TwoSlopeNorm(vmin=-500., vcenter=0, vmax=4000);
```
There are also various options for the ticks available in `matplotlib.ticker`.
"""
function set_colorbar(colorscale::ColorScale, v1, v2, data=[1.0f0])
   vmin, vmax = set_lim(Float32(v1), Float32(v2), data, colorscale)
   if colorscale == Linear
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      norm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
      ticks = matplotlib.ticker.LinearLocator(numticks=9)
   elseif colorscale == Log # logarithmic
      norm = matplotlib.colors.LogNorm(;vmin, vmax)
      ticks = matplotlib.ticker.LogLocator(base=10, subs=collect(0:9))
   else # symmetric log
      linthresh = 1.0
      logstep = 1

      logthresh = floor(Int, log10(linthresh))
      minlog = ceil(Int, log10(-vmin))
      maxlog = ceil(Int, log10(vmax))

      norm = matplotlib.colors.SymLogNorm(;linthresh, linscale=0.03, vmin, vmax, base=10)
      ticks = [ [-(10.0^x) for x in minlog:-logstep:logthresh]..., 0.0,
         [10.0^x for x in logthresh+1:logstep:maxlog]..., ]
   end

   norm, ticks
end

"Configure customized plot."
function set_plot(c, ax, pArgs::PlotArgs, ticks, addcolorbar)
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
      ax.spines[loc].set_linewidth(2.0)
   end

   ax.xaxis.set_tick_params(width=2.0, length=3)
   ax.yaxis.set_tick_params(width=2.0, length=3)
   return
end

"""
    vdfslice(meta::MetaVLSV, location::AbstractVector, ax=nothing; kwargs...)

Plot the 2D slice cut of phase space distribution function at `location` within velocity
range `limits`. If `ax===nothing`, plot on the current active axes.
# Optional arguments
- `unit::AxisUnit`: location unit in `SI`, `EARTH`.
- `unitv::String`: velocity unit in ("km/s", "m/s").
- `limits::Vector{Real}`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype`: symbol for choosing the slice type from `:xy`, `:xz`, `:yz`,
`:bperp`, `:bpar1`, `:bpar2`.
- `center`: symbol for setting the reference frame from `:bulk`, `:peak`.
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected!
- `vmin, vmax`: minimum and maximum VDF values for plotting.
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
- `flimit`: minimum VDF threshold for plotting.
- `kwargs...`: any valid keyword argument for hist2d.
"""
function vdfslice(meta::MetaVLSV, location, ax=nothing;
   limits=[-Inf32, Inf32, -Inf32, Inf32], verbose=false, species="proton",
   vmin=-Inf32, vmax=Inf32, unit=SI, unitv="km/s", slicetype=:default, vslicethick=0.0,
   center=:nothing, weight=:particle, flimit=-1.0,
   kwargs...)

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
    plotmesh(meta::MetaVLSV; projection="3d", origin=0.0, marker="+", kwargs...)

Plot mesh cell centers from axis view `projection`. `projection` should be either "3d", "x",
"y" or "z". `origin` is center of projection plane in the normal direction.
"""
function plotmesh(meta::MetaVLSV, ax=nothing; projection="3d", origin=0.0, marker="+",
   kwargs...)
   (;coordmin, coordmax, cellid) = meta
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
      ids = cellid
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
      s = ax.scatter(x1, x2, x3; marker, kwargs...)
   end
   s
end

"""
    pui(meta::MetaVLSV)
    pui(file::AbstractString)

Quick plotting via command line interactive selections.
"""
function pui(meta::MetaVLSV; suppress_output=false)

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