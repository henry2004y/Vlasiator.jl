# Plotting functionalities from Matplotlib.

using PyPlot

export plot, pcolormesh, pcolormeshslice, plot_vdf, streamplot, quiver, plotmesh

@static if matplotlib.__version__ >= "3.3"
   matplotlib.rc("image", cmap="turbo") # set default colormap
end

"Plotting arguments."
struct PlotArgs
   "data array size"
   sizes::Vector{Int}
   "plotting data range"
   plotrange::Vector{Float32}
   "cell IDs in the cut plane"
   idlist::Vector{Int}
   "mapping from original cell order to cut plane"
   indexlist::Vector{Int}
   "scale of data"
   colorscale::ColorScale
   "minimum data value"
   vmin::Float32
   "maximum data value"
   vmax::Float32
   "title"
   str_title::String
   "xlabel"
   strx::String
   "ylabel"
   stry::String
   "colorbar title"
   cb_title_use::String
end

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
   comp="xy", axisunit=RE, stride::Integer=10, kwargs...)

   X, Y, v1, v2 = set_vector(meta, var, comp, axisunit)

   if isnothing(ax) ax = plt.gca() end

   Xq  = @view X[1:stride:end, 1:stride:end]
   Yq  = @view Y[1:stride:end, 1:stride:end]
   v1q = @view v1[1:stride:end, 1:stride:end]
   v2q = @view v2[1:stride:end, 1:stride:end]

   ax.quiver(Xq, Yq, v1q, v2q; kwargs...)
end

function set_vector(meta::MetaVLSV, var, comp, axisunit)
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
   axisunit=RE, colorscale=Linear, addcolorbar=true, vmin=-Inf, vmax=Inf, kwargs...)

   if ndims(meta) == 3 || meta.maxamr > 0
      # check if origin and normal exist in kwargs
      normal = haskey(kwargs, :normal) ? kwargs.data.normal : :y
      origin = haskey(kwargs, :origin) ? kwargs.data.origin : 0.0
      kwargs = Base.structdiff(kwargs.data, (normal = normal, origin = origin))
      c = pcolormeshslice(meta, var, ax; op, axisunit, colorscale, addcolorbar, vmin, vmax,
         normal, origin, kwargs...)
      return c
   end

   pArgs = set_args(meta, var, axisunit, colorscale; vmin, vmax)

   x, y, data = plot_prep2d(meta, var, pArgs, op, axisunit)

   if var in ("fg_b", "fg_e", "vg_b_vol", "vg_e_vol") || endswith(var, "vg_v")
      rho_ = findfirst(endswith("rho"), meta.variable)
      if !isnothing(rho_)
         rho = readvariable(meta, meta.variable[rho_])
         rho = reshape(rho, pArgs.sizes[1], pArgs.sizes[2])
         mask = findall(==(0.0), rho')

         if ndims(data) == 2
            @inbounds data[mask] .= NaN
         else
            ind = CartesianIndices((pArgs.sizes[1], pArgs.sizes[2]))
            for m in mask
               @inbounds data[:, ind[m][1], ind[m][2]] .= NaN
            end
         end
      end
   end

   cnorm, cticks = set_colorbar(pArgs, data)

   if isnothing(ax) ax = plt.gca() end

   if colorscale != SymLog
      c = ax.pcolormesh(x, y, data; norm=cnorm, shading="auto", kwargs...)
   else
      c = ax.pcolormesh(x, y, data; norm=cnorm, cmap=matplotlib.cm.RdBu, shading="auto",
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

   pArgs = set_args(meta, var, axisunit, colorscale; normal, origin, vmin, vmax)

   @unpack sizes, plotrange, idlist, indexlist = pArgs

   data = readvariable(meta, var)

   if startswith(var, "fg_") # field quantities, fsgrid
      throw(ArgumentError("FS grid variable $var plotting in cut currently not supported!"))
   else # moments, dccrg grid
      # vlasov grid, AMR
      if ndims(data) == 1
         data = data[indexlist] # find required cells
      elseif ndims(data) == 2
         data = data[:,indexlist] # find required cells
      end

      # Create the plotting grid
      if ndims(data) == 1
         data = refineslice(meta, idlist, data, normal)
      elseif ndims(data) == 2
         if op in (:x, :y, :z, :1, :2, :3)
            if op in (:x, :1)
               slice = @view data[1,:]
            elseif op in (:y, :2)
               slice = @view data[2,:]
            elseif op in (:z, :3)
               slice = @view data[3,:]
            end
            data = refineslice(meta, idlist, slice, normal)
         elseif op == :mag
            datax = @views refineslice(meta, idlist, data[1,:], normal)
            datay = @views refineslice(meta, idlist, data[2,:], normal)
            dataz = @views refineslice(meta, idlist, data[3,:], normal)
            data = hypot.(datax, datay, dataz)
         end

      elseif ndims(data) == 3
         @error "not implemented yet!"
      end
   end

   x, y = get_axis(axisunit, plotrange, sizes)

   cnorm, cticks = set_colorbar(pArgs, data)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data'; norm=cnorm, shading="auto", kwargs...)

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   c
end

"Generate axis and data for 2D plotting."
function plot_prep2d(meta::MetaVLSV, var, pArgs::PlotArgs, op, axisunit::AxisUnit)
   @unpack sizes, plotrange = pArgs

   dataRaw = Vlasiator.getdata2d(meta, var)

   if ndims(dataRaw) == 3
      if op in (:x, :1)
         data = @view dataRaw[1,:,:]
      elseif op in (:y, :2)
         data = @view dataRaw[2,:,:]
      elseif op in (:z, :3)
         data = @view dataRaw[3,:,:]
      elseif op == :mag
         data = @views hypot.(dataRaw[1,:,:], dataRaw[2,:,:], dataRaw[3,:,:])
      end
   else
      data = dataRaw
   end

   x, y = Vlasiator.get_axis(axisunit, plotrange, sizes)

   x, y, data'
end

"Set plot-related arguments."
function set_args(meta::MetaVLSV, var, axisunit::AxisUnit, colorscale::ColorScale;
   normal::Symbol=:none, origin=0.0, vmin=-Inf, vmax=Inf)
   @unpack ncells, coordmin, coordmax = meta

   if normal == :x
      seq = @SVector [2,3]
      dir = 1
   elseif normal == :y || (ncells[2] == 1 && ncells[3] != 1) # polar
      seq = @SVector [1,3]
      dir = 2
   elseif normal == :z || (ncells[3] == 1 && ncells[2] != 1) # ecliptic
      seq = @SVector [1,2]
      dir = 3
   else
      throw(ArgumentError("1D data detected. Please use 1D plot functions."))
   end

   sizes = ncells[seq]
   plotrange = [coordmin[seq[1]], coordmax[seq[1]], coordmin[seq[2]], coordmax[seq[2]]]
   axislabels = ['X', 'Y', 'Z'][seq]

   if normal == :none
      idlist, indexlist = Int[], Int[]
   else
      idlist, indexlist = let sliceoffset = abs(coordmin[dir]) + origin
         getslicecell(meta, sliceoffset, dir, coordmin[dir], coordmax[dir])
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^meta.maxamr # data needs to be refined later

   unitstr = axisunit == RE ? L"$R_E$" : L"$m$"
   strx = axislabels[1]*"["*unitstr*"]"
   stry = axislabels[2]*"["*unitstr*"]"

   str_title = @sprintf "t= %4.1fs" meta.time

   datainfo = readvariablemeta(meta, var)

   cb_title_use = !isempty(datainfo.variableLaTeX) ?
      datainfo.variableLaTeX * " ["*datainfo.unitLaTeX*"]" : ""

   PlotArgs(sizes, plotrange, idlist, indexlist, colorscale,
      vmin, vmax, str_title, strx, stry, cb_title_use)
end

"Set colorbar norm and ticks."
function set_colorbar(pArgs::PlotArgs, data=[1.0])
   @unpack colorscale, vmin, vmax = pArgs
   if colorscale == Linear
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(v1, v2)
      cnorm = matplotlib.colors.BoundaryNorm(levels, ncolors=256, clip=true)
      ticks = matplotlib.ticker.LinearLocator(numticks=9)
   elseif colorscale == Log # logarithmic
      datapositive = data[data .> 0.0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax

      cnorm = matplotlib.colors.LogNorm(vmin=v1, vmax=v2)
      ticks = matplotlib.ticker.LogLocator(base=10, subs=collect(0:9))
   else # symmetric log
      linthresh = 1.0
      logstep = 1

      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax

      logthresh = floor(Int, log10(linthresh))
      minlog = ceil(Int, log10(-v1))
      maxlog = ceil(Int, log10(v2))

      cnorm = matplotlib.colors.SymLogNorm(;linthresh, linscale=0.03, vmin=v1, vmax=v2,
         base=10)
      ticks = [ [-(10.0^x) for x in minlog:-logstep:logthresh]..., 0.0,
         [10.0^x for x in logthresh+1:logstep:maxlog]..., ]
   end

   cnorm, ticks
end


"Configure customized plot."
function set_plot(c, ax, pArgs::PlotArgs, cticks, addcolorbar)
   @unpack str_title, strx, stry, cb_title_use = pArgs

   if addcolorbar
      cb = colorbar(c; ax, ticks=cticks, fraction=0.046, pad=0.04)
      !isempty(cb_title_use) && cb.ax.set_ylabel(cb_title_use, fontsize=14)
      cb.outline.set_linewidth(1.0)
   end

   ax.set_title(str_title, fontsize=14, fontweight="bold")
   ax.set_xlabel(strx, fontsize=14)
   ax.set_ylabel(stry, fontsize=14)
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
   vxsize = vmesh.vblocks[1] * vmesh.vblock_size[1]
   vysize = vmesh.vblocks[2] * vmesh.vblock_size[2]

   vxmin, vxmax = vmesh.vmin[1], vmesh.vmax[1]
   vymin, vymax = vmesh.vmin[2], vmesh.vmax[2]

   cellsize = (vxmax - vxmin) / vxsize # this assumes cubic vspace grid!

   unit == RE && (location .*= Vlasiator.Re)

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

   # Set normal direction
   if slicetype == :nothing
      if ncells[2] == 1 && ncells[3] == 1 # 1D, select xz
         slicetype = :xz
      elseif ncells[2] == 1 # polar
         slicetype = :xz
      elseif ncells[3] == 1 # ecliptic
         slicetype == :xy
      end
   elseif slicetype in (:bperp, :bpar, :bpar1)
      # If necessary, find magnetic field
      if hasvariable(meta, "B_vol")
         B = readvariable(meta, "B_vol", cidNearest)
      elseif hasvariable(meta, "vg_b_vol")
         B = readvariable(meta, "vg_b_vol", cidNearest)
      end
      BxV = B × Vbulk
      if slicetype == :bperp # slice in b_perp1/b_perp2
         sliceNormal = B ./ norm(B)
         strx = L"$v_{B \times V}$ "
         stry = L"$v_{B \times (B \times V)}$ "
      elseif slicetype == :bpar1 # slice in b_parallel/b_perp1 plane
         sliceNormal = B × BxV
         sliceNormal ./= norm(sliceNormal)
         strx = L"$v_{B}$ "
         stry = L"$v_{B \times V}$ "
      else # slice in b_parallel/b_perp2 plane
         sliceNormal = BxV ./ norm(BxV)
         strx = L"$v_{B}$ "
         stry = L"$v_{B \times (B \times V)}$ "
      end
   end

   if slicetype == :xy
      v1 = V[1,:]
      v2 = V[2,:]
      vnormal = V[3,:]
      sliceNormal = [0., 0., 1.]
      strx = "vx [km/s]"
      stry = "vy [km/s]"
   elseif slicetype == :yz
      v1 = V[2,:]
      v2 = V[3,:]
      vnormal = V[1,:]
      strx = "vy [km/s]"
      stry = "vz [km/s]"
   elseif slicetype == :xz
      v1 = V[1,:]
      v2 = V[3,:]
      vnormal = V[2,:]
      sliceNormal = [0., 1., 0.]
      strx = "vx [km/s]"
      stry = "vz [km/s]"
   elseif slicetype ∈ (:Bperp, :Bpar, :Bpar1)
      #TODO: NOT working yet!
      if slicetype == :Bperp
         v1 = Vrot2[1,:] # the X axis of the slice is BcrossV=perp1
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[3,:] # the Z axis of the slice is B
      elseif slicetype == :Bpar
         v1 = Vrot2[3,:] # the X axis of the slice is B
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[1,:] # the Z axis of the slice is -BcrossV=perp1
      elseif slicetype == :Bpara1
         v1 = Vrot2[3,:] # the X axis of the slice is B
         v2 = Vrot2[1,:] # the Y axis of the slice is BcrossV=perp1
         vnormal = Vrot2[2,:] # the Z axis of the slice is Bcross(BcrossV)=perp2
      end
   end

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
         (vxmin < v1 < vxmax) & (vymin < v2 < vymax)
   else
      ind_ = @. (vxmin < v1 < vxmax) & (vymin < v2 < vymax)
   end

   # [m/s] --> [km/s]
   unitfactor = 1e3
   v1, v2, fw = v1[ind_]./unitfactor, v2[ind_]./unitfactor, fw[ind_]

   isinf(fmin) && (fmin = minimum(fw))
   isinf(fmax) && (fmax = maximum(fw))

   verbose && @info "Active f range is $fmin, $fmax"

   if isnothing(ax) ax = plt.gca() end

   cnorm = matplotlib.colors.LogNorm(vmin=fmin, vmax=fmax)

   rx = LinRange(vxmin/unitfactor, vxmax/unitfactor, vxsize+1)
   ry = LinRange(vymin/unitfactor, vymax/unitfactor, vysize+1)

   h = ax.hist2d(v1, v2, bins=(rx, ry), weights=fw, norm=cnorm)

   ax.set_title(str_title, fontsize=14, fontweight="bold")
   ax.set_xlabel(strx, fontsize=14, weight="black")
   ax.set_ylabel(stry, fontsize=14, weight="black")
   ax.set_aspect("equal")
   ax.grid(color="grey", linestyle="-")

   cb = colorbar(h[4]; ax=ax, fraction=0.046, pad=0.04)
   cb_title = cb.ax.set_ylabel("f(v)", fontsize=14)

   if slicetype in ("bperp", "bpar", "bpar1")
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