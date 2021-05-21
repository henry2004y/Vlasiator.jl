# Plotting functionalities from Matplotlib.

using PyPlot, Printf, LaTeXStrings
using LinearAlgebra: norm, ×

import PyPlot: plot, quiver, streamplot, pcolormesh
export plot, pcolormesh, pcolormeshslice, plot_vdf, streamplot, quiver, plotmesh

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
   "maximum refinement level"
   maxreflevel::Int8
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
   "colormap used in colorbar"
   cmap::ColorMap
   "colorbar title"
   cb_title_use::String
end

"""
    plot(meta, var, ax=nothing; kwargs)

Plot `var` from `meta` of 1D VLSV data. If `ax===nothing`, plot on the current active axes.
"""
function plot(meta::MetaData, var, ax=nothing; kwargs...)
   if hasvariable(meta, var)
      data = readvariable(meta, var)
   else
      data = Vlasiator.variables_predefined[var](meta)
   end

   x = LinRange(meta.xmin, meta.xmax, meta.xcells)

   if isnothing(ax) ax = plt.gca() end

   ax.plot(x, data; kwargs...)
end

"""
    streamplot(meta, var, ax=nothing; comp="xy", axisunit=RE, kwargs...)

Wrapper over Matplotlib's streamplot function. The `comp` option can take a subset of "xyz"
in any order. `axisunit` can be chosen from `RE, SI`.
The keyword arguments can be any valid Matplotlib arguments into streamplot.
"""
function streamplot(meta::MetaData, var, ax=nothing; comp="xy", axisunit=RE, kwargs...)

   if occursin("x", comp)
      v1_ = 1
      if occursin("y", comp)
         v2_ = 2
         sizes = [meta.xcells, meta.ycells]
         plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
      else
         v2_ = 3
         sizes = [meta.xcells, meta.zcells]
         plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      end
   else
      v1_, v2_ = 2, 3
      sizes = [meta.ycells, meta.zcells]
      plotrange = [meta.ymin, meta.ymax, meta.zmin, meta.zmax]
   end

   if hasvariable(meta, var)
      data = readvariable(meta, var)
   else
      data = Vlasiator.variables_predefined[var](meta)
   end

   if startswith(var, "fg_")
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   else # vlasov grid
      @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required!"
      data = reshape(data, 3, sizes...)
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   end

   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   # Be careful about array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   if isnothing(ax) ax = plt.gca() end

   streamplot(X, Y, v1, v2; kwargs...)
end

"""
    quiver(meta, var, ax=nothing; comp="xy", axisunit=RE, kwargs...)

Wrapper over Matplotlib's quiver function. If `ax===nothing`, plot on the current active
axes. The `comp` option can take a subset of "xyz" in any order. `axisunit` can be chosen
from `RE, SI`.
The keyword arguments can be any valid Matplotlib arguments into quiver.
"""
function quiver(meta::MetaData, var, ax=nothing; comp="xy", axisunit=RE, kwargs...)

   if occursin("x", comp)
      v1_ = 1
      if occursin("y", comp)
         v2_ = 2
         sizes = [meta.xcells, meta.ycells]
         plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
      else
         v2_ = 3
         sizes = [meta.xcells, meta.zcells]
         plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      end
   else
      v1_, v2_ = 2, 3
      sizes = [meta.ycells, meta.zcells]
      plotrange = [meta.ymin, meta.ymax, meta.zmin, meta.zmax]
   end

   if hasvariable(meta, var)
      data = readvariable(meta, var)
   else
      data = Vlasiator.variables_predefined[var](meta)
   end

   if startswith(var, "fg_")
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   else # vlasov grid
      @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required!"
      data = reshape(data, 3, sizes...)
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   end

   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   # Be careful about array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   if isnothing(ax) ax = plt.gca() end

   ax.quiver(X, Y, v1, v2; kwargs...)
end

"""
    pcolormesh(meta::MetaData, var, ax=nothing;
       op=:mag, axisunit=RE, colorscale=Log, vmin=-Inf, vmax=Inf, addcolorbar=true)

Plot a variable using pseudocolor from 2D VLSV data.
If `ax` is provided, then it will plot on that axes.
If 3D or AMR grid detected, it will pass arguments to [`pcolormeshslice`](@ref).

# Optional arguments
- `op::Symbol`: the component of a vector to plot, chosen from `:mag, :x, :y, :z`.
- `axisunit::AxisUnit`: the unit of axis ∈ `RE, SI`.
- `colorscale::ColorScale`: whether to use linear scale for data.
- `vmin::Float`: minimum data range. Set to maximum of data if not specified. 
- `vmax::Float`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`pcolormesh(meta, var)`

`pcolormesh(meta, var, axisunit=SI)`

`pcolormesh(data, func, colorscale=Linear)`
"""
function pcolormesh(meta::MetaData, var, ax=nothing;
   op=:mag, axisunit=RE, colorscale=Log, addcolorbar=true, vmin=-Inf, vmax=Inf, kwargs...)

   if showdimension(meta) == 3 || getmaxamr(meta) > 0
      # check if origin and normal exist in kwargs
      normal = haskey(kwargs, :normal) ? kwargs.data.normal : :y
      origin = haskey(kwargs, :origin) ? kwargs.data.origin : 0.0
      kwargs = Base.structdiff(kwargs.data, (normal = normal, origin = origin))
      c = pcolormeshslice(meta, var, ax; op, axisunit, colorscale, addcolorbar, vmin, vmax,
         normal, origin, kwargs...)
      return c
   end

   pArgs = set_args(meta, var, axisunit, colorscale; normal=:none, vmin, vmax)

   x, y, data = plot_prep2d(meta, var, pArgs, op, axisunit)

   cnorm, cticks = set_colorbar(pArgs, data)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data; norm=cnorm, cmap=pArgs.cmap, shading="auto", kwargs...)

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   return c
end

"""
    pcolormeshslice(meta, var, ax=nothing; kwargs...)

Plot pseudocolor var on a 2D slice of 3D vlsv data.
If `ax` is provided, then it will plot on that axes.
It would be easier to call [`pcolormesh`](@ref).

# Optional arguments
- `op::Symbol`: the component of a vector to plot, chosen from `:mag, :x, :y, :z`.
- `origin::Float`: center of slice plane in the normal direction.
- `normal::Symbol`: the normal direction of cut plane, chosen from `:x, :y, :z`.
- `axisunit::AxisUnit`: the unit of axis ∈ `RE, SI`.
- `colorscale::ColorScale`: color scale for data ∈ (`Linear`, `Log`)
- `vmin::Real`: minimum data range. Set to maximum of data if not specified. 
- `vmax::Real`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`pcolormeshslice(meta, var)`

`pcolormeshslice(meta, var, op=:z, origin=1.0, normal=:x)`

`pcolormeshslice(data, func, colorscale=Log)`
"""
function pcolormeshslice(meta::MetaData, var, ax=nothing; op::Symbol=:mag, origin=0.0,
   normal::Symbol=:y, axisunit::AxisUnit=RE, colorscale::ColorScale=Log, addcolorbar=true,
   vmin::Real=-Inf, vmax::Real=Inf, kwargs...)

   pArgs = set_args(meta, var, axisunit, colorscale; normal, origin, vmin, vmax)

   maxreflevel = pArgs.maxreflevel
   sizes = pArgs.sizes
   plotrange = pArgs.plotrange
   idlist, indexlist = pArgs.idlist, pArgs.indexlist

   if hasvariable(meta, var)
      data = readvariable(meta, var)
   else
      data = Vlasiator.variables_predefined[var](meta)
   end

   if startswith(var, "fg_") # field quantities, fsgrid

   else # moments, dccrg grid
      # vlasov grid, AMR
      if ndims(data) == 1
         data = data[indexlist] # find required cells
      elseif ndims(data) == 2
         data = data[:,indexlist] # find required cells
      end

      # Create the plotting grid
      if ndims(data) == 1
         data = refineslice(meta, idlist, data, maxreflevel, normal)
      elseif ndims(data) == 2
         if op in (:x, :y, :z)
            if op == :x
               slice = @view data[1,:]
            elseif op == :y
               slice = @view data[2,:]
            elseif op == :z
               slice = @view data[3,:]
            end
            data = refineslice(meta, idlist, slice, maxreflevel, normal)
         elseif op == :mag
            datax = @views refineslice(meta, idlist, data[1,:], maxreflevel, normal)
            datay = @views refineslice(meta, idlist, data[2,:], maxreflevel, normal)
            dataz = @views refineslice(meta, idlist, data[3,:], maxreflevel, normal)
            data = hypot.(datax, datay, dataz)
         end

      elseif ndims(data) == 3
         @error "not implemented yet!"
      end
   end

   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   cnorm, cticks = set_colorbar(pArgs, data)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data'; norm=cnorm, cmap=pArgs.cmap, shading="auto", kwargs...)

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   c
end

"Generate axis and data for 2D plotting."
function plot_prep2d(meta, var, pArgs, op, axisunit::AxisUnit)

   sizes, plotrange = pArgs.sizes, pArgs.plotrange

   if hasvariable(meta, var)
      dataRaw = readvariable(meta, var)
   else
      dataRaw = Vlasiator.variables_predefined[var](meta)
   end

   if ndims(dataRaw) == 1 || (ndims(dataRaw) == 2 && size(dataRaw)[1] == 1)
      data = reshape(dataRaw, sizes[1], sizes[2])
   else
      if ndims(dataRaw) == 2
         dataRaw = reshape(dataRaw, 3, sizes...)
      end
      if op == :x
         data = dataRaw[1,:,:]
      elseif op == :y
         data = dataRaw[2,:,:]
      elseif op == :z
         data = dataRaw[3,:,:]
      elseif op == :mag
         data = hypot.(dataRaw[1,:,:], dataRaw[2,:,:], dataRaw[3,:,:])
      end
   end

   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   x, y, data'
end

"Set plot-related arguments."
function set_args(meta, var, axisunit::AxisUnit, colorscale::ColorScale;
   normal::Symbol=:z, origin=0.0, vmin=-Inf, vmax=Inf)

   maxreflevel = getmaxamr(meta)

   if normal == :x
      sizes = [meta.ycells, meta.zcells]
      plotrange = [meta.ymin, meta.ymax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.xmin) + origin
      axislabels = ['Y', 'Z']

      idlist, indexlist = getslicecell(meta, sliceoffset, maxreflevel;
         xmin=meta.xmin, xmax=meta.xmax)
   elseif normal == :y
      sizes = [meta.xcells, meta.zcells]
      plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.ymin) + origin
      axislabels = ['X', 'Z']

      idlist, indexlist = getslicecell(meta, sliceoffset, maxreflevel;
         ymin=meta.ymin, ymax=meta.ymax)
   elseif normal == :z
      sizes = [meta.xcells, meta.ycells]
      plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
      sliceoffset = abs(meta.zmin) + origin
      axislabels = ['X', 'Y']

      idlist, indexlist = getslicecell(meta, sliceoffset, maxreflevel;
         zmin=meta.zmin, zmax=meta.zmax)
   else
      idlist = Int64[]
      indexlist = Int64[]

      if meta.ycells == 1 && meta.zcells != 1 # polar
         plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
         sizes = [meta.xcells, meta.zcells]
         PLANE = "XZ"
         axislabels = ['X', 'Z']
      elseif meta.zcells == 1 && meta.ycells != 1 # ecliptic
         plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
         sizes = [meta.xcells, meta.ycells]
         PLANE = "XY"
         axislabels = ['X', 'Y']
      else # 1D
         @error "1D data detected. Please use 1D plot functions."
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^maxreflevel # data needs to be refined later

   unitstr = axisunit == RE ? "R_E" : "m"
   strx = latexstring(axislabels[1]*"["*unitstr*"]")
   stry = latexstring(axislabels[2]*"["*unitstr*"]")

   if hasparameter(meta, "t")
      timesim = readparameter(meta, "t")
      str_title = @sprintf "t= %4.1fs" timesim
   elseif hasparameter(meta, "time")
      timesim = readparameter(meta, "time")
      str_title = @sprintf "t= %4.1fs" timesim
   else
      str_title = ""
   end

   cmap = matplotlib.cm.turbo

   datainfo = readvariablemeta(meta, var)

   cb_title_use = datainfo.variableLaTeX
   cb_title_use *= ",["*datainfo.unitLaTeX*"]"

   PlotArgs(sizes, plotrange, idlist, indexlist, maxreflevel, colorscale,
      vmin, vmax, str_title, strx, stry, cmap, cb_title_use)
end

"Set colorbar norm and ticks."
function set_colorbar(pArgs, data)

   if pArgs.colorscale == Log # Logarithmic plot
      datapositive = data[data .> 0.0]
      if isempty(datapositive)
         throw(DomainError(data, "Nonpositive data detected: use linear scale instead!"))
      end
      vmin = isinf(pArgs.vmin) ? minimum(datapositive) : pArgs.vmin
      vmax = isinf(pArgs.vmax) ? maximum(data) : pArgs.vmax

      cnorm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
      ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   elseif pArgs.colorscale == Linear
      vmin = isinf(pArgs.vmin) ? minimum(data) : pArgs.vmin
      vmax = isinf(pArgs.vmax) ? maximum(data) : pArgs.vmax
      nticks = 7
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      cnorm = matplotlib.colors.BoundaryNorm(levels, ncolors=pArgs.cmap.N,
         clip=true)
      ticks = range(vmin, vmax, length=nticks)
   end

   cnorm, ticks
end


"Configure customized plot."
function set_plot(c, ax, pArgs, cticks, addcolorbar)

   str_title, strx, stry, cb_title_use = pArgs.str_title, pArgs.strx, pArgs.stry,
      pArgs.cb_title_use

   if addcolorbar
      cb = colorbar(c; ax, ticks=cticks, fraction=0.046, pad=0.04)
      cb_title = cb.ax.set_ylabel(cb_title_use, fontsize=14)
      cb.outline.set_linewidth(1.0)
   end

   ax.set_title(str_title, fontsize=14, fontweight="bold")
   ax.set_xlabel(strx, fontsize=14, weight="black")
   ax.set_ylabel(stry, fontsize=14, weight="black")
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
- `slicetype`: string for choosing the slice type from "xy", "xz", "yz", "bperp", "bpar",
"bpar1".
- `center`: string for setting the reference frame from "bulk", "peak".
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected! 
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
"""
function plot_vdf(meta, location, ax=nothing; limits=[-Inf, Inf, -Inf, Inf],
   verbose=false, pop="proton", fmin=-Inf, fmax=Inf, unit::AxisUnit=SI, slicetype="xy",
   vslicethick=0.0, center="", weight::Symbol=:particle, fThreshold=-1.0)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells

   vmesh = meta.meshes[pop]
   vxsize = vmesh.vxblocks * vmesh.vxblock_size
   vysize = vmesh.vyblocks * vmesh.vyblock_size
   vzsize = vmesh.vzblocks * vmesh.vzblock_size
   vxmin, vxmax = vmesh.vxmin, vmesh.vxmax
   vymin, vymax = vmesh.vymin, vmesh.vymax
   vzmin, vzmax = vmesh.vzmin, vmesh.vzmax
   cellsize = (vxmax - vxmin) / vxsize # this assumes cubic vspace grid!

   unit == RE && (location ./= Vlasiator.Re)

   if pop == "proton"
      if !Vlasiator.hasname(meta.footer, "BLOCKIDS", "proton")
         if Vlasiator.hasname(meta.footer, "BLOCKIDS", "avgs") # old versions
            pop = "avgs"
         else
            @error "Unable to detect population "*pop
         end
      end
   elseif !Vlasiator.hasname(meta.footer, "BLOCKIDS", pop)
      @error "Unable to detect population "*pop
   end

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

   for f in ("fsaved", "vg_f_saved")
      if hasvariable(meta, f) &&
         readvariable(meta, f, cidNearest) != 1.0
         @error "VDF not found in the given cell!"
      end
   end
      
   vcellids, vcellf = readvcells(meta, cidNearest; pop)

   V = getvcellcoordinates(meta, vcellids; pop)

   if center == "bulk" # center with bulk velocity
      verbose && @info "Transforming to plasma frame"
      V -= Vbulk
   elseif center == "peak" # center on highest f-value
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

   if hasparameter(meta, "t")
      timesim = readparameter(meta, "t")
      str_title = @sprintf "t= %4.1fs" timesim
   elseif hasparameter(meta, "time")
      timesim = readparameter(meta, "time")
      str_title = @sprintf "t= %4.1fs" timesim
   else
      str_title = ""
   end

   # Set normal direction
   if ysize == 1 && zsize == 1 # 1D, select xz
      slicetype = "xz"
      sliceNormal = [0., 1., 0.]
      strx = "vx [km/s]"
      stry = "vz [km/s]"
   elseif ysize == 1 && slicetype == "xz" # polar
      sliceNormal = [0., 1., 0.]
      strx = "vx [km/s]"
      stry = "vz [km/s]"
   elseif zsize == 1 && slicetype == "xy" # ecliptic
      sliceNormal = [0., 0., 1.]
      strx = "vx [km/s]"
      stry = "vy [km/s]"
   elseif slicetype in ("bperp", "bpar", "bpar1")
      # If necessary, find magnetic field
      if hasvariable(meta, "B_vol")
         B = readvariable(meta, "B_vol", cidNearest)
      elseif hasvariable(meta, "vg_b_vol")
         B = readvariable(meta, "vg_b_vol", cidNearest)
      end
      BxV = B × Vbulk
      if slicetype == "bperp" # slice in b_perp1/b_perp2
         sliceNormal = B ./ norm(B)
         strx = L"$v_{B \times V}$ "
         stry = L"$v_{B \times (B \times V)}$ "
      elseif slicetype == "bpar1" # slice in b_parallel/b_perp1 plane
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

   if slicetype == "xy"
      v1 = V[1,:]
      v2 = V[2,:]
      vnormal = V[3,:]
   elseif slicetype == "yz"
      v1 = V[2,:]
      v2 = V[3,:]
      vnormal = V[1,:]
   elseif slicetype == "xz"
      v1 = V[1,:]
      v2 = V[3,:]
      vnormal = V[2,:]
   elseif slicetype ∈ ("Bperp", "Bpar", "Bpar1")
      #hyzhou: NOT working yet!
      if slicetype == "Bperp"
         v1 = Vrot2[1,:] # the X axis of the slice is BcrossV=perp1
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[3,:] # the Z axis of the slice is B
      elseif slicetype == "Bpar"
         v1 = Vrot2[3,:] # the X axis of the slice is B
         v2 = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[1,:] # the Z axis of the slice is -BcrossV=perp1
      elseif slicetype == "Bpara1"
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
   cmap = matplotlib.cm.turbo

   rx = LinRange(vxmin/unitfactor, vxmax/unitfactor, vxsize+1)
   ry = LinRange(vymin/unitfactor, vymax/unitfactor, vysize+1)

   h = ax.hist2d(v1, v2, bins=(rx, ry), weights=fw, norm=cnorm, cmap=cmap)

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
function plotmesh(meta::MetaData, ax=nothing; projection="3d", origin=0.0, marker="+",
   kwargs...)

   maxreflevel = getmaxamr(meta)

   if projection == "x"
      sliceoffset = abs(meta.xmin) + origin

      ids, _ = getslicecell(meta, sliceoffset, maxreflevel;
         xmin=meta.xmin, xmax=meta.xmax)
   elseif projection == "y"
      sliceoffset = abs(meta.ymin) + origin

      ids, _ = getslicecell(meta, sliceoffset, maxreflevel;
         ymin=meta.ymin, ymax=meta.ymax)
   elseif projection == "z"
      sliceoffset = abs(meta.zmin) + origin

      ids, _ = getslicecell(meta, sliceoffset, maxreflevel;
         zmin=meta.zmin, zmax=meta.zmax)
   else
      ids = meta.cellid
   end

   centers = Matrix{Float32}(undef, 3, length(ids))
   for (i, id) in enumerate(ids)
      centers[:,i] = getcellcoordinates(meta, id)
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