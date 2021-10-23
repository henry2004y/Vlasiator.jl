# Plotting functionalities from Matplotlib.

using PyPlot

export plot, pcolormesh, pcolormeshslice, plot_vdf, streamplot, quiver, plotmesh

@static if matplotlib.__version__ >= "3.3"
   matplotlib.rc("image", cmap="turbo") # set default colormap
end

matplotlib.rc("font", size=14)
matplotlib.rc("xtick", labelsize=10) 
matplotlib.rc("ytick", labelsize=10) 

"""
    plot(meta, var, ax=nothing; kwargs)

Plot `var` from `meta` of 1D VLSV data. If `ax===nothing`, plot on the current active axes.
"""
function PyPlot.plot(meta::MetaVLSV, var, ax=nothing; kwargs...)
   data = readvariable(meta, var)

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   if isnothing(ax) ax = plt.gca() end

   ax.plot(x, data; kwargs...)
end

"""
    streamplot(meta, var, ax=nothing; comp="xy", axisunit=RE, kwargs...)

Wrapper over Matplotlib's streamplot function. The `comp` option can take a subset of "xyz"
in any order. `axisunit` can be chosen from `RE, SI`.
The keyword arguments can be any valid Matplotlib arguments into streamplot.

# Optional arguments
- `comp`: a subset of "xyz" in any order.
- `axisunit`: chosen from RE and SI.
"""
function PyPlot.streamplot(meta::MetaVLSV, var::AbstractString, ax=nothing;
   comp="xy", axisunit=RE, kwargs...)

   X, Y, v1, v2 = set_vector(meta, var, comp, axisunit)

   if isnothing(ax) ax = plt.gca() end

   streamplot(X, Y, v1, v2; kwargs...)
end

"""
    quiver(meta, var, ax=nothing; comp="xy", axisunit=RE, stride=10, kwargs...)

Wrapper over Matplotlib's quiver function. If `ax===nothing`, plot on the current active
axes. The `comp` option can take a subset of "xyz" in any order. `axisunit` can be chosen
from `RE, SI`.
The keyword arguments can be any valid Matplotlib arguments into quiver.

# Optional arguments
- `comp`: a subset of "xyz" in any order.
- `axisunit`: chosen from RE and SI.
- `stride::Integer`: arrow strides in number of cells.
"""
function PyPlot.quiver(meta::MetaVLSV, var::AbstractString, ax=nothing;
   comp="xy", axisunit::AxisUnit=RE, stride::Integer=10, kwargs...)

   X, Y, v1, v2 = set_vector(meta, var, comp, axisunit)

   if isnothing(ax) ax = plt.gca() end

   Xq  = @view X[1:stride:end, 1:stride:end]
   Yq  = @view Y[1:stride:end, 1:stride:end]
   v1q = @view v1[1:stride:end, 1:stride:end]
   v2q = @view v2[1:stride:end, 1:stride:end]

   ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

function set_vector(meta::MetaVLSV, var, comp, axisunit::AxisUnit)
   @unpack ncells, coordmin, coordmax = meta
   if occursin("x", comp)
      v1_ = 1
      if occursin("y", comp)
         v2_ = 2
         sizes = [ncells[1], ncells[2]]
         plotrange = [coordmin[1], coordmax[1], coordmin[2], coordmax[2]]
      else
         v2_ = 3
         sizes = [ncells[1], ncells[3]]
         plotrange = [coordmin[1], coordmax[1], coordmin[3], coordmax[3]]
      end
   else
      v1_, v2_ = 2, 3
      sizes = [ncells[2], ncells[3]]
      plotrange = [coordmin[2], coordmax[2], coordmin[3], coordmax[3]]
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
       op=:mag, axisunit=RE, colorscale=Linear, vmin=-Inf, vmax=Inf, addcolorbar=true,
       kwargs...)

Plot a variable using pseudocolor from 2D VLSV data.
If `ax` is provided, then it will plot on that axes.
If 3D or AMR grid detected, it will pass arguments to [`pcolormeshslice`](@ref).

# Optional arguments
- `op::Symbol`: the component of a vector, chosen from `:mag, :x, :y, :z, :1, :2, :3`.
- `axisunit::AxisUnit`: the unit of axis ∈ `RE, SI`.
- `colorscale::ColorScale`: Linear, Log, or SymLog.
- `vmin::Float`: minimum data range. Set to maximum of data if not specified.
- `vmax::Float`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`pcolormesh(meta, var)`

`pcolormesh(meta, var, axisunit=SI)`

`pcolormesh(data, func, colorscale=Log)`
"""
function PyPlot.pcolormesh(meta::MetaVLSV, var::AbstractString, ax=nothing; op=:mag,
   axisunit::AxisUnit=RE, colorscale::ColorScale=Linear, addcolorbar=true,
   vmin=-Inf, vmax=Inf, kwargs...)

   if ndims(meta) == 3 || meta.maxamr > 0
      # check if origin and normal exist in kwargs
      normal = haskey(kwargs, :normal) ? kwargs.data.normal : :y
      origin = haskey(kwargs, :origin) ? kwargs.data.origin : 0.0
      kwargs = Base.structdiff(kwargs.data, (normal = normal, origin = origin))
      c = pcolormeshslice(meta, var, ax; op, axisunit, colorscale, addcolorbar, vmin, vmax,
         normal, origin, kwargs...)
      return c
   end

   pArgs = set_args(meta, var, axisunit)

   x, y = get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
   data = prep2d(meta, var, op)'

   if var in ("fg_b", "fg_e", "vg_b_vol", "vg_e_vol") || endswith(var, "vg_v")
      rho_ = findfirst(endswith("rho"), meta.variable)
      if !isnothing(rho_)
         rho = readvariable(meta, meta.variable[rho_])
         rho = reshape(rho, pArgs.sizes[1], pArgs.sizes[2])
         mask = findall(==(0.0), rho')

         if ndims(data) == 2
            @inbounds data[mask] .= NaN
         else
            ind = CartesianIndices((pArgs.sizes[2], pArgs.sizes[1]))
            for m in mask
               @inbounds data[:, ind[m][1], ind[m][2]] .= NaN
            end
         end
      end
   end

   cnorm, cticks = set_colorbar(colorscale, vmin, vmax, data)

   if isnothing(ax) ax = plt.gca() end

   if colorscale != SymLog
      c = ax.pcolormesh(x, y, data; norm=cnorm, shading="nearest", kwargs...)
   else
      c = ax.pcolormesh(x, y, data; norm=cnorm, cmap=matplotlib.cm.RdBu, shading="nearest",
         kwargs...)
   end

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   return c
end

"""
    pcolormeshslice(meta, var, ax=nothing; kwargs...)

Plot pseudocolor var on a 2D slice of 3D vlsv data.
If `ax` is provided, then it will plot on that axes.
It would be easier to call [`pcolormesh`](@ref), since it auto-detects dimension.

# Optional arguments
- `op::Symbol`: the component of a vector, chosen from `:mag, :x, :y, :z, :1, :2, :3`.
- `origin::Float`: center of slice plane in the normal direction.
- `normal::Symbol`: the normal direction of cut plane, chosen from `:x, :y, :z`.
- `axisunit::AxisUnit`: the unit of axis ∈ `RE, SI`.
- `colorscale::ColorScale`: color scale for data ∈ (`Linear`, `Log`, `SymLog`)
- `vmin::Real`: minimum data range. Set to maximum of data if not specified.
- `vmax::Real`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`pcolormeshslice(meta, var)`

`pcolormeshslice(meta, var, op=:z, origin=1.0, normal=:x)`

`pcolormeshslice(data, func, colorscale=Log)`
"""
function pcolormeshslice(meta::MetaVLSV, var::AbstractString, ax=nothing; op::Symbol=:mag,
   origin=0.0, normal::Symbol=:y, axisunit::AxisUnit=RE, colorscale::ColorScale=Linear,
   addcolorbar=true, vmin::Real=-Inf, vmax::Real=Inf, kwargs...)

   pArgs = set_args(meta, var, axisunit; normal, origin)

   data = prep2dslice(meta, var, normal, op, pArgs)'
   x, y = get_axis(axisunit, pArgs.plotrange, pArgs.sizes)

   cnorm, cticks = set_colorbar(colorscale, vmin, vmax, data)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data; norm=cnorm, shading="nearest", kwargs...)

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   c
end

"Set colorbar norm and ticks in a given range `v1` to `v2` for `data` in `colorscale`."
function set_colorbar(colorscale::ColorScale, v1, v2, data=[1.0])
   vmin, vmax = set_lim(v1, v2, data, colorscale)
   if colorscale == Linear
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      cnorm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
      ticks = matplotlib.ticker.LinearLocator(numticks=9)
   elseif colorscale == Log # logarithmic
      cnorm = matplotlib.colors.LogNorm(;vmin, vmax)
      ticks = matplotlib.ticker.LogLocator(base=10, subs=collect(0:9))
   else # symmetric log
      linthresh = 1.0
      logstep = 1

      logthresh = floor(Int, log10(linthresh))
      minlog = ceil(Int, log10(-vmin))
      maxlog = ceil(Int, log10(vmax))

      cnorm = matplotlib.colors.SymLogNorm(;linthresh, linscale=0.03, vmin, vmax, base=10)
      ticks = [ [-(10.0^x) for x in minlog:-logstep:logthresh]..., 0.0,
         [10.0^x for x in logthresh+1:logstep:maxlog]..., ]
   end

   cnorm, ticks
end

"Configure customized plot."
function set_plot(c, ax, pArgs::PlotArgs, cticks, addcolorbar)
   @unpack str_title, strx, stry, cb_title = pArgs

   if addcolorbar
      cb = colorbar(c; ax, ticks=cticks, fraction=0.04, pad=0.02)
      !isempty(cb_title) && cb.ax.set_ylabel(cb_title)
      cb.ax.tick_params(direction="in")
   end

   ax.set_title(str_title, fontweight="bold")
   ax.set_xlabel(strx)
   ax.set_ylabel(stry)
   ax.set_aspect("equal")

   # Set border line widths
   for loc in ("left", "bottom", "right", "top")
      edge = get(ax.spines, loc, nothing)
      edge.set_linewidth(2.0)
   end

   ax.xaxis.set_tick_params(width=2.0, length=3)
   ax.yaxis.set_tick_params(width=2.0, length=3)
   return
end

"""
    plot_vdf(meta, location, ax=nothing; kwargs...)

Plot the 2D slice cut of phase space distribution function at `location` within velocity
range `limits`. If `ax===nothing`, plot on the current active axes.
# Optional arguments
- `unit::AxisUnit`: axis unit in `SI`, `RE`.
- `limits::Vector{Real}`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype`: symbol for choosing the slice type from :xy, :xz, :yz, :bperp, :bpar, :bpar1.
- `center`: symbol for setting the reference frame from :bulk, :peak.
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected!
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
- `kwargs...`: any valid keyword argument for hist2d.
"""
function plot_vdf(meta::MetaVLSV, location, ax=nothing; limits=[-Inf, Inf, -Inf, Inf],
   verbose=false, pop="proton", fmin=-Inf, fmax=Inf, unit::AxisUnit=SI,
   slicetype::Symbol=:nothing, vslicethick=0.0, center::Symbol=:nothing,
   weight::Symbol=:particle, fThreshold=-1.0, kwargs...)

   @unpack ncells = meta
   if haskey(meta.meshes, pop)
      vmesh = meta.meshes[pop]
   else
      throw(ArgumentError("Unable to detect population $pop"))
   end

   unit == RE && (location .*= Re)

   # Calculate cell ID from given coordinates
   cidReq = getcell(meta, location)
   cidNearest = getnearestcellwithvdf(meta, cidReq)

   if verbose
      @info "Original coordinates : $location"
      @info "Original cell        : $(getcellcoordinates(meta, cidReq))"
      @info "Nearest cell with VDF: $(getcellcoordinates(meta, cidNearest))"
   end

   x, y, z = getcellcoordinates(meta, cidNearest)
   verbose && @info "cellid $cidNearest, x = $x, y = $y, z = $z"

   # Set normal direction
   if slicetype == :nothing
      if ncells[2] == 1 && ncells[3] == 1 # 1D, select xz
         slicetype = :xz
      elseif ncells[2] == 1 # polar
         slicetype = :xz
      elseif ncells[3] == 1 # ecliptic
         slicetype = :xy
      end
   end

   if slicetype == :xy
      dir1, dir2, dir3 = 1, 2, 3
      sliceNormal = [0., 0., 1.]
   elseif slicetype == :xz
      dir1, dir2, dir3 = 1, 3, 2
      sliceNormal = [0., 1., 0.]
   elseif slicetype == :yz
      dir1, dir2, dir3 = 2, 3, 1
      sliceNormal = [1., 0., 0.]
   elseif slicetype in (:bperp, :bpar, :bpar1)
      if hasvariable(meta, "B_vol")
         B = readvariable(meta, "B_vol", cidNearest)
      elseif hasvariable(meta, "vg_b_vol")
         B = readvariable(meta, "vg_b_vol", cidNearest)
      end
      BxV = B × Vbulk
      if slicetype == :bperp # slice in b_perp1/b_perp2
         sliceNormal = B ./ norm(B)
      elseif slicetype == :bpar1 # slice in b_parallel/b_perp1 plane
         sliceNormal = B × BxV
         sliceNormal ./= norm(sliceNormal)
      else # slice in b_parallel/b_perp2 plane
         sliceNormal = BxV ./ norm(BxV)
      end
   end
   
   v1size = vmesh.vblocks[dir1] * vmesh.vblock_size[dir1]
   v2size = vmesh.vblocks[dir2] * vmesh.vblock_size[dir2]

   v1min, v1max = vmesh.vmin[dir1], vmesh.vmax[dir1]
   v2min, v2max = vmesh.vmin[dir2], vmesh.vmax[dir2]

   @assert (v1max - v1min) / v1size ≈ (v2max - v2min) / v2size "Noncubic vgrid detected!"
   cellsize = (v1max - v1min) / v1size

   # Extracts Vbulk
   if hasvariable(meta, "moments")
      # This should be a restart file
      Vbulk = readvariable(meta, "restart_V", cidNearest)
   elseif hasvariable(meta, pop*"/vg_v")
      # multipop v5 bulk file
      Vbulk = readvariable(meta, pop*"/vg_v", cidNearest)
   elseif hasvariable(meta, pop*"/V")
      # multipop bulk file
      Vbulk = readvariable(meta, pop*"/V", cidNearest)
   else
      # regular bulk file, currently analysator supports pre- and
      # post-multipop files with "V"
      Vbulk = readvariable(meta, "V", cidNearest)
   end

   vcellids, vcellf = readvcells(meta, cidNearest; pop)

   V = getvcellcoordinates(meta, vcellids; pop)

   if center == :bulk # centered with bulk velocity
      verbose && @info "Transforming to plasma frame"
      V -= Vbulk
   elseif center == :peak # centered on highest f-value
      peakindex = argmax(vcellf)
      Vpeak = V[:,peakindex]
      V -= Vpeak
      verbose && "Plot in frame of peak f-value, travelling at speed $Vpeak"
   end

   # Set sparsity threshold
   if hasvariable(meta, pop*"/EffectiveSparsityThreshold")
      fThreshold = readvariable(meta,
         pop*"/EffectiveSparsityThreshold", cidNearest)
   elseif hasvariable(meta, pop*"/vg_effectivesparsitythreshold")
      fThreshold = readvariable(meta,
         pop+"/vg_effectivesparsitythreshold", cidNearest)
   else
      verbose && @info "Using a default f threshold value of 1e-16."
      fThreshold = 1e-16
   end

   # Drop all velocity cells which are below the sparsity threshold
   fselect_ = vcellf .≥ fThreshold
   f = vcellf[fselect_]
   V = V[:,fselect_]

   str_title = @sprintf "t= %4.1fs" meta.time

   if slicetype ∈ (:xy, :yz, :xz)
      v1 = V[dir1,:]
      v2 = V[dir2,:]
      vnormal = V[dir3,:]
      strx, stry = getindex(["vx", "vy", "vz"], [dir1, dir2])
   elseif slicetype ∈ (:Bperp, :Bpar, :Bpar1)
      #TODO: NOT working yet!
      if slicetype == :Bperp
         v1 = Vrot2[1,:] # the X axis of the slice is BcrossV=perp1
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[3,:] # the Z axis of the slice is B
         strx = L"$v_{B \times V}$"
         stry = L"$v_{B \times (B \times V)}$"
      elseif slicetype == :Bpar
         v1 = Vrot2[3,:] # the X axis of the slice is B
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[1,:] # the Z axis of the slice is -BcrossV=perp1
         strx = L"$v_{B}$"
         stry = L"$v_{B \times V}$"
      elseif slicetype == :Bpar1
         v1 = Vrot2[3,:] # the X axis of the slice is B
         v2 = Vrot2[1,:] # the Y axis of the slice is BcrossV=perp1
         vnormal = Vrot2[2,:] # the Z axis of the slice is Bcross(BcrossV)=perp2
         strx = L"$v_{B}$"
         stry = L"$v_{B \times (B \times V)}$"
      end
   end

   strx *= " [km/s]"
   stry *= " [km/s]"

   # Weights using particle flux or phase-space density
   fw = weight == :flux ? f*norm([v1, v2, vnormal]) : f

   if verbose
      if vslicethick > 0
         @info "Performing slice with a counting thickness of $vslicethick"
      else
         @info "Projecting total VDF to a single plane"
      end
   end

   if vslicethick < 0 # Trying to set a proper value automatically
      # Assure that the slice cut through at least 1 velocity cell
      if any(sliceNormal .== 1.0)
         vslicethick = cellsize
      else # Assume cubic vspace grid, add extra space
         vslicethick = cellsize*(√3+0.05)
      end
   end

   # Select cells which are within slice area
   if vslicethick > 0.0
      ind_ = @. (abs(vnormal) ≤ 0.5*vslicethick) &
         (v1min < v1 < v1max) & (v2min < v2 < v2max)
   else
      ind_ = @. (v1min < v1 < v1max) & (v2min < v2 < v2max)
   end

   # [m/s] --> [km/s]
   unitvfactor = 1e3
   v1, v2, fw = v1[ind_]./unitvfactor, v2[ind_]./unitvfactor, fw[ind_]

   isinf(fmin) && (fmin = minimum(fw))
   isinf(fmax) && (fmax = maximum(fw))

   verbose && @info "Active f range is $fmin, $fmax"

   if isnothing(ax) ax = plt.gca() end

   cnorm = matplotlib.colors.LogNorm(vmin=fmin, vmax=fmax)

   r1 = LinRange(v1min/unitvfactor, v1max/unitvfactor, v1size+1)
   r2 = LinRange(v2min/unitvfactor, v2max/unitvfactor, v2size+1)

   h = ax.hist2d(v1, v2, bins=(r1, r2), weights=fw, norm=cnorm)

   ax.set_title(str_title, fontweight="bold")
   ax.set_xlabel(strx, weight="black")
   ax.set_ylabel(stry, weight="black")
   ax.set_aspect("equal")
   ax.grid(color="grey", linestyle="-")

   cb = colorbar(h[4]; ax=ax, fraction=0.046, pad=0.02)
   cb_title = cb.ax.set_ylabel("f(v)")

   if slicetype in (:bperp, :bpar, :bpar1)
      # Draw vector of magnetic field direction
   end
   plt.tight_layout()

   h[4] # h[1] is 2D data, h[2] is x axis, h[3] is y axis
end

"""
    plotmesh(meta; projection="3d", origin=0.0, marker="+", kwargs...)

Plot mesh cell centers from axis view `projection`. `projection` should be either "3d", "x",
"y" or "z". `origin` is center of projection plane in the normal direction.
"""
function plotmesh(meta::MetaVLSV, ax=nothing; projection="3d", origin=0.0, marker="+",
   kwargs...)
   @unpack coordmin, coordmax, cellid = meta
   if projection == "x"
      sliceoffset = origin - coordmin[1]
      ids, _ = getslicecell(meta, sliceoffset, 1, coordmin[1], coordmax[1])
   elseif projection == "y"
      sliceoffset = origin - coordmin[2]
      ids, _ = getslicecell(meta, sliceoffset, 2, coordmin[2], coordmax[2])
   elseif projection == "z"
      sliceoffset = origin - coordmin[3]
      ids, _ = getslicecell(meta, sliceoffset, 3, coordmin[3], coordmax[3])
   else
      ids = cellid
   end

   centers = Matrix{Float32}(undef, 3, length(ids))
   for (i, id) in enumerate(ids)
      @inbounds centers[:,i] = getcellcoordinates(meta, id)
   end

   if isnothing(ax) ax = plt.gca() end

   if projection == "x"
      s = ax.scatter(centers[2,:], centers[3,:]; marker, kwargs...)
   elseif projection == "y"
      s = ax.scatter(centers[1,:], centers[3,:]; marker, kwargs...)
   elseif projection == "z"
      s = ax.scatter(centers[1,:], centers[2,:]; marker, kwargs...)
   else
      s = ax.scatter(centers[1,:], centers[2,:], centers[3,:]; marker, kwargs...)
   end
   s
end