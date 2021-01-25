# Vlasiator plotting in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu 12/03/2020

using PyPlot, Printf, LaTeXStrings

export streamplot, plot_pcolormesh, plot_colormap3dslice

import PyPlot:streamplot

"""
    streamplot(meta::MetaData, var; comp="xy", axisunit="Re", kwargs...)

Wrapper over Matplotlib's stream line function. The `comp` option can take a
subset of "xyz" in any order.
"""
function streamplot(meta, var; comp="xy", axisunit="Re", kwargs...)

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

   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required to plot streamlines!"

   if startswith(var, "fg_") # fsgrid
      
   else # vlasov grid
      data = reshape(data[:,meta.cellIndex], 3, sizes...)
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   end

   x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Re
   y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Re

   # Be careful about array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   c = streamplot(X, Y, v1, v2; kwargs...)
end


"""
    plot_pcolormesh(meta::MetaData, var; op="mag", axisunit="Re", islinear=false)

Plot a 2D pseudocolor var from vlsv.

`plot_pcolormesh(meta, var)`

`plot_pcolormesh(meta, var, axisunit="SI")`

`plot_pcolormesh(data, func, islinear=false)`
"""
function plot_pcolormesh(meta, var; op="mag", axisunit="Re", islinear=false)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin 
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   dx = meta.dx # cell size is equal in x,y,z for now

   # Check if ecliptic or polar run
   if ysize == 1 && zsize != 1
      plotrange = [xmin, xmax, zmin, zmax]
      sizes = [xsize, zsize]
      PLANE = "XZ"
      axislabels = ['X', 'Z']
   elseif zsize == 1 && ysize != 1
      plotrange = [xmin, xmax, ymin, ymax]
      sizes = [xsize, ysize]
      PLANE = "XY"
      axislabels = ['X', 'Y']
   elseif ysize == 1 && zsize == 1

   end

   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   if startswith(var, "fg_") # fsgrid
      
   else # vlasov grid
      if ndims(data) == 1 || (ndims(data) == 2 && size(data)[1] == 1)       
         data = reshape(data[meta.cellIndex], sizes[1], sizes[2])
      elseif ndims(data) == 2
         data = reshape(data[:,meta.cellIndex], 3, sizes...)
         if op == "x"
            data = data[1,:,:]
         elseif op == "y"
            data = data[2,:,:]
         elseif op == "z"
            data = data[3,:,:]
         elseif startswith("mag",op)
            @. data = sqrt(data[1,:,:] ^2 + data[2,:,:]^2 + data[3,:,:]^2)
         end
      elseif ndims(data) == 3

      else
         @error "Error in reshaping data $(var)!"
      end
   end

   x, y, str_title, strx, stry, cmap, norm, ticks, cb_title_use = 
      set_args(meta, var, axisunit, islinear, axislabels, plotrange, sizes, data)

   fig, ax = subplots()

   c = ax.pcolormesh(x, y, data', norm=norm, cmap=cmap, shading="auto")

   set_plot(fig, ax, c, str_title, strx, stry, ticks, cb_title_use)

   return c
end

"""
    plot_colormap3dslice(meta::MetaData, var (...))

Plot pseudocolor var on a 2D slice of 3D vlsv data.

`plot_colormap3dslice(meta, var)`

`plot_colormap3dslice(meta, var, op="z", origin=1.0, normal="x")`

`plot_colormap3dslice(data, func, islinear=false)`
"""
function plot_colormap3dslice(meta, var; op="mag", origin=0.0, normal="y",
   axisunit="Re", islinear=false)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin 
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax

   maxreflevel = get_max_amr_level(meta)

   if normal == "x"
      normal_D = [1,0,0]
      sizes = [ysize, zsize]
      plotrange = [ymin, ymax, zmin, zmax]
      sliceoffset = abs(xmin) + origin
      axislabels = ['Y','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         xmin=xmin, xmax=xmax)
   elseif normal == "y"
      normal_D = [0,1,0]
      sizes = [xsize, zsize]
      plotrange = [xmin, xmax, zmin, zmax]
      sliceoffset = abs(ymin) + origin
      axislabels = ['X','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         ymin=ymin, ymax=ymax)
   elseif normal == "z"
      normal_D = [0,0,1]
      sizes = [xsize, ysize]
      plotrange = [xmin, xmax, ymin, ymax]
      sliceoffset = abs(zmin) + origin
      axislabels = ['X','Y']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         zmin=zmin, zmax=zmax)
   end

   # Scale the sizes to the heighest refinement level
   sizes *= 2^maxreflevel

   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   if startswith(var, "fg_") # field quantities, fsgrid

   else # moments, dccrg grid
      # vlasov grid, AMR
      if ndims(data) == 1
         data = data[meta.cellIndex] # sorted
         data = data[indexlist] # find required cells
      elseif ndims(data) == 2
         data = data[:,meta.cellIndex] # sorted
         data = data[:,indexlist] # find required cells
      end

      # Create the plotting grid
      if ndims(data) == 1
         data = refine_data(meta, idlist, data, maxreflevel, normal)
      elseif ndims(data) == 2
         if op == "x"
            data = refine_data(meta, idlist, data[1,:], maxreflevel, normal)
         elseif op == "y"
            data = refine_data(meta, idlist, data[2,:], maxreflevel, normal)
         elseif op == "z"
            data = refine_data(meta, idlist, data[3,:], maxreflevel, normal)
         elseif startswith("mag",op)
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

   x, y, str_title, strx, stry, cmap, norm, ticks, cb_title_use = 
      set_args(meta, var, axisunit, islinear, axislabels, plotrange, sizes, data)

   fig, ax = subplots()

   c = ax.pcolormesh(x, y, data', norm=norm, cmap=cmap, shading="auto")

   set_plot(fig, ax, c, str_title, strx, stry, ticks, cb_title_use)

   return c
end

"Set plot-related arguments."
function set_args(meta, var, axisunit, islinear, axislabels, plotrange, sizes, data)

   if axisunit == "Re"
      x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Re
      y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Re      
   else
      x = range(plotrange[1], plotrange[2], length=sizes[1])
      y = range(plotrange[3], plotrange[4], length=sizes[2])
   end

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

   if !islinear
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
      norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=true)
      ticks = range(vmin, vmax, length=nticks)
   end

   datainfo = read_variable_info(meta, var)

   cb_title_use = datainfo.variableLaTeX
   data_unit = datainfo.unitLaTeX
   cb_title_use *= ",["*data_unit*"]"

   return x, y, str_title, strx, stry, cmap, norm, ticks, cb_title_use
end

"Draw customized plot."
function set_plot(fig, ax, c, str_title, strx, stry, ticks, cb_title_use)

   cb = fig.colorbar(c, ticks=ticks)

   ax.set_title(str_title, fontsize=14, fontweight="bold")
   ax.set_xlabel(strx, fontsize=14, weight="black")
   ax.set_ylabel(stry, fontsize=14, weight="black")
   ax.set_aspect("equal")

   for axis in ["top","bottom","left","right"]
      ax.spines[axis].set_linewidth(2.0)
   end
   ax.xaxis.set_tick_params(width=2.0,length=3)
   ax.yaxis.set_tick_params(width=2.0,length=3)
   
   cb_title = cb.ax.set_title(cb_title_use, fontsize=14, fontweight="bold")
   cb.outline.set_linewidth(2.0)
end