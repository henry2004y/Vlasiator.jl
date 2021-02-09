# Vlasiator plotting in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Printf, LaTeXStrings

"Plotting arguments."
struct PlotArgs
   sizes::Vector{Int}
   plotrange::Vector{Float32}
   idlist::Vector{Int}
   indexlist::Vector{Int}
   maxreflevel::Int8
   islinear::Bool
   str_title::String
   strx::String
   stry::String
   cmap::ColorMap
   cb_title_use::String
end

"""
    streamline(meta::MetaData, var; comp="xy", axisunit="Re", kwargs...)

Wrapper over Matplotlib's streamplot function. The `comp` option can take a
subset of "xyz" in any order.
"""
function streamline(meta, var; comp="xy", axisunit="Re", kwargs...)

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

   if var in keys(Vlasiator.variables_predefined)
      data = Vlasiator.variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required to plot streamlines!"

   if startswith(var, "fg_") # fsgrid
      
   else # vlasov grid
      data = reshape(data, 3, sizes...)
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   end

   x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Vlasiator.Re
   y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Vlasiator.Re

   # Be careful about array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   c = streamplot(X, Y, v1, v2; kwargs...)
end

"""
    plot_pcolormesh(meta::MetaData, var, ax=nothing; op=:mag, axisunit="Re",
       islinear=false)

Plot a variable using pseudocolor from 2D VLSV data. If `ax` is provided, then
it tries to plot into that axes.

`plot_pcolormesh(meta, var)`

`plot_pcolormesh(meta, var, axisunit="SI")`

`plot_pcolormesh(data, func, islinear=false)`
"""
function plot_pcolormesh(meta, var, ax=nothing; op=:mag, axisunit="Re",
   islinear=false, addcolorbar=true)

   pArgs = set_args(meta, var, axisunit, islinear)

   x, y, data = plot_prep2d(meta, var, pArgs, op, axisunit)

   norm, cticks = set_colorbar(data, pArgs)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data, norm=norm, cmap=pArgs.cmap, shading="auto")

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   return c
end

"""
    plot_colormap3dslice(meta::MetaData, var, ax=nothing (...))

Plot pseudocolor var on a 2D slice of 3D vlsv data. If `ax` is provided, then
it tries to plot into that axes.

`plot_colormap3dslice(meta, var)`

`plot_colormap3dslice(meta, var, op=:z, origin=1.0, normal=:x)`

`plot_colormap3dslice(data, func, islinear=false)`
"""
function plot_colormap3dslice(meta, var, ax=nothing; op=:mag, origin=0.0,
   normal=:y, axisunit="Re", islinear=false, addcolorbar=true)

   pArgs = set_args(meta, var, axisunit, islinear; normal, origin)

   maxreflevel = pArgs.maxreflevel
   sizes = pArgs.sizes
   plotrange = pArgs.plotrange
   idlist, indexlist = pArgs.idlist, pArgs.indexlist

   if var in keys(Vlasiator.variables_predefined)
      data = Vlasiator.variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
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
         data = refine_data(meta, idlist, data, maxreflevel, normal)
      elseif ndims(data) == 2
         if op == :x
            data = refine_data(meta, idlist, data[1,:], maxreflevel, normal)
         elseif op == :y
            data = refine_data(meta, idlist, data[2,:], maxreflevel, normal)
         elseif op == :z
            data = refine_data(meta, idlist, data[3,:], maxreflevel, normal)
         elseif op == :mag
            datax = refine_data(meta, idlist, data[1,:], maxreflevel, normal)
            datay = refine_data(meta, idlist, data[2,:], maxreflevel, normal)
            dataz = refine_data(meta, idlist, data[3,:], maxreflevel, normal)
            data = @. sqrt(datax^2 + datay^2 + dataz^2)
         end

      elseif ndims(data) == 3
         @error "not implemented yet!"
      else
         @error "Dimension error in constructing 2D AMR slice!"
      end
   end

   if axisunit == "Re"
      x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Vlasiator.Re
      y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Vlasiator.Re      
   else
      x = range(plotrange[1], plotrange[2], length=sizes[1])
      y = range(plotrange[3], plotrange[4], length=sizes[2])
   end

   norm, cticks = set_colorbar(data, pArgs)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data', norm=norm, cmap=pArgs.cmap, shading="auto")

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   return c
end

"Generate axis and data for 2D plotting."
function plot_prep2d(meta, var, pArgs, op, axisunit)

   sizes, plotrange = pArgs.sizes, pArgs.plotrange

   if var in keys(Vlasiator.variables_predefined)
      data = Vlasiator.variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   if ndims(data) == 1 || (ndims(data) == 2 && size(data)[1] == 1)       
      data = reshape(data, sizes[1], sizes[2])
   else
      if ndims(data) == 2
         data = reshape(data, 3, sizes...)
      end
      if op == :x
         data = data[1,:,:]
      elseif op == :y
         data = data[2,:,:]
      elseif op == :z
         data = data[3,:,:]
      elseif op == :mag
         data = @. sqrt(data[1,:,:] ^2 + data[2,:,:]^2 + data[3,:,:]^2)
      end
   end

   if axisunit == "Re"
      x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Vlasiator.Re
      y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Vlasiator.Re
   else
      x = range(plotrange[1], plotrange[2], length=sizes[1])
      y = range(plotrange[3], plotrange[4], length=sizes[2])
   end

   return x, y, data'
end

"Set plot-related arguments."
function set_args(meta, var, axisunit, islinear; normal=:y, origin=0.0)

   maxreflevel = get_max_amr_level(meta)

   if normal == :x
      normal_D = [1,0,0]
      sizes = [meta.ycells, meta.zcells]
      plotrange = [meta.ymin, meta.ymax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.xmin) + origin
      axislabels = ['Y','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         xmin=meta.xmin, xmax=meta.xmax)
   elseif normal == :y
      normal_D = [0,1,0]
      sizes = [meta.xcells, meta.zcells]
      plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.ymin) + origin
      axislabels = ['X','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         ymin=meta.ymin, ymax=meta.ymax)
   elseif normal == :z
      normal_D = [0,0,1]
      sizes = [meta.xcells, meta.ycells]
      plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
      sliceoffset = abs(meta.zmin) + origin
      axislabels = ['X','Y']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         zmin=meta.zmin, zmax=meta.zmax)
   else
      idlist = Int64[]
      indexlist = Int64[]
      # Check if ecliptic or polar run
      if meta.ycells == 1 && meta.zcells != 1
         plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
         sizes = [meta.xcells, meta.zcells]
         PLANE = "XZ"
         axislabels = ['X', 'Z']
      elseif meta.zcells == 1 && meta.ycells != 1
         plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
         sizes = [meta.xcells, meta.ycells]
         PLANE = "XY"
         axislabels = ['X', 'Y']
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^maxreflevel # data needs to be refined later (WIP)

   strx = latexstring(axislabels[1]*"["*axisunit*"]")
   stry = latexstring(axislabels[2]*"["*axisunit*"]")

   if has_parameter(meta, "t")
      timesim = read_parameter(meta, "t")
      str_title = @sprintf "t= %4.1fs" timesim
   elseif has_parameter(meta, "time")
      timesim = read_parameter(meta, "time")
      str_title = @sprintf "t= %4.1fs" timesim
   else
      str_title = ""
   end

   cmap = matplotlib.cm.turbo

   datainfo = read_variable_info(meta, var)

   cb_title_use = datainfo.variableLaTeX
   cb_title_use *= ",["*datainfo.unitLaTeX*"]"

   PlotArgs(sizes, plotrange, idlist, indexlist, maxreflevel, islinear,
      str_title, strx, stry, cmap, cb_title_use)
end

function set_colorbar(data, pArgs)

   if !pArgs.islinear
      # Logarithmic plot
      vmin = minimum(data[data .> 0.0])
      vmax = maximum(data)
      
      norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
      ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   else
      vmin = minimum(data)
      vmax = maximum(data)
      nticks = 7
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      norm = matplotlib.colors.BoundaryNorm(levels, ncolors=pArgs.cmap.N, clip=true)
      ticks = range(vmin, vmax, length=nticks)
   end

   return norm, ticks
end


"Configure customized plot."
function set_plot(c, ax, pArgs, cticks, addcolorbar)

   str_title, strx, stry, cb_title_use = pArgs.str_title,
      pArgs.strx, pArgs.stry, pArgs.cb_title_use

   if addcolorbar
      cb = colorbar(c, ax=ax, ticks=cticks, fraction=0.046, pad=0.04)
      cb_title = cb.ax.set_ylabel(cb_title_use, fontsize=14)
      cb.outline.set_linewidth(1.0)
   end

   ax.set_title(str_title, fontsize=14, fontweight="bold")
   ax.set_xlabel(strx, fontsize=14, weight="black")
   ax.set_ylabel(stry, fontsize=14, weight="black")
   ax.set_aspect("equal")

   for axis in ["top","bottom","left","right"]
      ax.spines[axis].set_linewidth(2.0)
   end
   ax.xaxis.set_tick_params(width=2.0,length=3)
   ax.yaxis.set_tick_params(width=2.0,length=3)
end