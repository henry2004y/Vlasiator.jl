# Vlasiator plotting in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu

using Vlasiator, PyPlot, Printf, LaTeXStrings
import LinearAlgebra: norm, ×

export plot_pcolormesh, plot_colormap3dslice, plot_vdf, streamline

"Plotting arguments."
struct PlotArgs
   sizes::Vector{Int}
   plotrange::Vector{Float32}
   idlist::Vector{Int}        # cell IDs in the cut plane
   indexlist::Vector{Int}     # mapping from original cell order to cut plane
   maxreflevel::Int8          # maximum refinement level
   islinear::Bool             # linear scale data
   vmin::Float32              # minimum data value 
   vmax::Float32              # maximum data value
   str_title::String
   strx::String
   stry::String
   cmap::ColorMap
   cb_title_use::String
end

"""
    streamline(meta::MetaData, var; comp="xy", axisunit="Re", kwargs...)

Wrapper over Matplotlib's streamplot function. The `comp` option can take a
subset of "xyz" in any order. `axisunit` can be chosen from `"Re", "SI"`.
The keyword arguments can be any valid Matplotlib arguments into streamplot.
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

   if startswith(var, "fg_")
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   else # vlasov grid
      @assert ndims(data) == 2 && size(data,1) == 3 "Vector data required!"
      data = reshape(data, 3, sizes...)
      v1 = data[v1_,:,:]'
      v2 = data[v2_,:,:]'
   end

   if axisunit == "Re"
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   # Be careful about array ordering difference between Julia and Python!
   X = [i for _ in y, i in x]
   Y = [j for j in y, _ in x]

   c = streamplot(X, Y, v1, v2; kwargs...)
end

"""
    plot_pcolormesh(meta::MetaData, var, ax=nothing; op=:mag, axisunit="Re",
       islinear=false, vmin=-Inf, vmax=Inf, addcolorbar=true)

Plot a variable using pseudocolor from 2D VLSV data. If `ax` is provided, then
it will plot on that axes.

# Optional arguments
- `op::Symbol`: the component of a vector to plot, chosen from `:mag, :x, :y, :z`.
- `axisunit::String`: the unit of axis, `"Re", "SI"`.
- `islinear::Bool`: whether to use linear scale for data.
- `vmin::Float`: minimum data range. Set to maximum of data if not specified. 
- `vmax::Float`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`plot_pcolormesh(meta, var)`

`plot_pcolormesh(meta, var, axisunit="SI")`

`plot_pcolormesh(data, func, islinear=false)`
"""
function plot_pcolormesh(meta, var, ax=nothing; op=:mag, axisunit="Re",
   islinear=false, addcolorbar=true, vmin=-Inf, vmax=Inf)

   pArgs = set_args(meta, var, axisunit, islinear; normal=:none, vmin, vmax)

   x, y, data = plot_prep2d(meta, var, pArgs, op, axisunit)

   cnorm, cticks = set_colorbar(pArgs, data)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data, norm=cnorm, cmap=pArgs.cmap, shading="auto")

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   c
end

"""
    plot_colormap3dslice(meta::MetaData, var, ax=nothing (...))

Plot pseudocolor var on a 2D slice of 3D vlsv data. If `ax` is provided, then
it will plot on that axes.

# Optional arguments
- `op::Symbol`: the component of a vector to plot, chosen from `:mag, :x, :y, :z`.
- `origin::Float`: center of slice plane in the normal direction.
- `normal::Symbol`: the normal direction of cut plane, chosen from `:x, :y, :z`.
- `axisunit::String`: the unit of axis, `"Re", "SI"`.
- `islinear::Bool`: whether to use linear scale for data.
- `vmin::Float`: minimum data range. Set to maximum of data if not specified. 
- `vmax::Float`: maximum data range. Set to minimum of data if not specified.
- `addcolorbar::Bool`: whether to add a colorbar to the colormesh.

`plot_colormap3dslice(meta, var)`

`plot_colormap3dslice(meta, var, op=:z, origin=1.0, normal=:x)`

`plot_colormap3dslice(data, func, islinear=false)`
"""
function plot_colormap3dslice(meta, var, ax=nothing; op=:mag, origin=0.0,
   normal=:y, axisunit="Re", islinear=false, addcolorbar=true,
   vmin=-Inf, vmax=Inf)

   pArgs = set_args(meta, var, axisunit, islinear;
      normal, origin, vmin=-Inf, vmax=Inf)

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
         if op in (:x, :y, :z)
            if op == :x
               slice = @view data[1,:]
            elseif op == :y
               slice = @view data[2,:]
            elseif op == :z
               slice = @view data[3,:]
            end
            data = refine_data(meta, idlist, slice, maxreflevel, normal)
         elseif op == :mag
            datax = @views refine_data(meta, idlist, data[1,:], maxreflevel, normal)
            datay = @views refine_data(meta, idlist, data[2,:], maxreflevel, normal)
            dataz = @views refine_data(meta, idlist, data[3,:], maxreflevel, normal)
            data = hypot.(datax, datay, dataz)
         end

      elseif ndims(data) == 3
         @error "not implemented yet!"
      else
         @error "Dimension error in constructing 2D AMR slice!"
      end
   end

   if axisunit == "Re"
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   cnorm, cticks = set_colorbar(data, pArgs)

   if isnothing(ax) ax = plt.gca() end

   c = ax.pcolormesh(x, y, data', norm=cnorm, cmap=pArgs.cmap, shading="auto")

   set_plot(c, ax, pArgs, cticks, addcolorbar)

   c
end

"Generate axis and data for 2D plotting."
function plot_prep2d(meta, var, pArgs, op, axisunit)

   sizes, plotrange = pArgs.sizes, pArgs.plotrange

   if var in keys(Vlasiator.variables_predefined)
      dataRaw = Vlasiator.variables_predefined[var](meta)
   else
      dataRaw = read_variable(meta, var)
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

   if axisunit == "Re"
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   x, y, data'
end

"Set plot-related arguments."
function set_args(meta, var, axisunit, islinear; normal=:z, origin=0.0,
   vmin=-Inf, vmax=Inf)

   maxreflevel = get_max_amr_level(meta)

   if normal == :x
      sizes = [meta.ycells, meta.zcells]
      plotrange = [meta.ymin, meta.ymax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.xmin) + origin
      axislabels = ['Y','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         xmin=meta.xmin, xmax=meta.xmax)
   elseif normal == :y
      sizes = [meta.xcells, meta.zcells]
      plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      sliceoffset = abs(meta.ymin) + origin
      axislabels = ['X','Z']

      idlist, indexlist = getSliceCellID(meta, sliceoffset, maxreflevel,
         ymin=meta.ymin, ymax=meta.ymax)
   elseif normal == :z
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
      vmin, vmax, str_title, strx, stry, cmap, cb_title_use)
end

"Return colorbar norm and ticks."
function set_colorbar(pArgs, data)

   if !pArgs.islinear # Logarithmic plot
      vmin = isinf(pArgs.vmin) ? minimum(data[data .> 0.0]) : pArgs.vmin
      vmax = isinf(pArgs.vmax) ? maximum(data) : pArgs.vmax

      cnorm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
      ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   else
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

"""
    plot_vdf(meta, location; kwargs...)

Plot the 2D slice cut of phase space distribution function at `location` within
velocity range `limits`.
# Optional arguments
- `limits`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype`: string for choosing the slice type from "xy", "xz",
"yz", "bperp", "bpar", "bpar1".
- `center`: string for setting the reference frame from "bulk", "peak".
- `vslicethick`: setting the velocity space slice thickness in the normal
direction. If set to 0, the whole distribution along the normal direction is
projected onto a plane. Currently this is only meaningful when `center` is set
such that a range near the bulk/peak normal velocity is selected! 
- `weight`: symbol for choosing distribution weights from phase space
density or particle flux.
"""
function plot_vdf(meta, location; limits=[-Inf, Inf, -Inf, Inf], ax=nothing,
   verbose=false, pop="proton", fmin=-Inf, fmax=Inf, unit="SI", slicetype="xy",
   vslicethick=0.0, center="0", weight=:particle, fThreshold=-1.0)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells

   vmesh = meta.meshes[pop]
   vxsize = vmesh.vxblocks * vmesh.vxblock_size
   vysize = vmesh.vyblocks * vmesh.vyblock_size
   vzsize = vmesh.vzblocks * vmesh.vzblock_size
   vxmin, vxmax = vmesh.vxmin, vmesh.vxmax
   vymin, vymax = vmesh.vymin, vmesh.vymax
   vzmin, vzmax = vmesh.vzmin, vmesh.vzmax
   cellsize = (vxmax - vxmin) / vxsize # this assumes cubic vspace grid!

   unit == "Re" && (location ./= Vlasiator.Re)

   if pop == "proton"
      if !Vlasiator.has_name(meta.footer, "BLOCKIDS", "proton")
         if Vlasiator.has_name(meta.footer, "BLOCKIDS", "avgs") # old versions
            pop = "avgs"
         else
            @error "Unable to detect population "*pop
         end
      end
   elseif !Vlasiator.has_name(meta.footer, "BLOCKIDS", pop)
      @error "Unable to detect population "*pop
   end

   # Calculate cell IDs from given coordinates        
   xReq = @view location[1,:]
   yReq = @view location[2,:]
   zReq = @view location[3,:]

   cellids = Int[]
   for i = 1:size(location, 2)
      cidReq = get_cellid(meta, [xReq[i], yReq[i], zReq[i]])
      cidNearest = getNearestCellWithVspace(meta, cidReq)

      if verbose
         @info "Point: $i out of $(size(location, 2)) requested"
         @info "Original coordinates : $(xReq[i]), $(yReq[i]), $(zReq[i])"
         @info "Original cell        : $(get_cell_coordinates(meta, cidReq))"
         @info "Nearest cell with VDF: $(get_cell_coordinates(meta, cidNearest))"
      end
      push!(cellids, cidNearest)
   end
   sort!(cellids); unique!(cellids)

   for cid in cellids
      x, y, z = get_cell_coordinates(meta, cid)
      verbose && @info "cellid $cid, x = $x, y = $y, z = $z"

      # Extracts Vbulk
      if has_variable(meta.footer, "moments")
         # This should be a restart file
         Vbulk = read_variable_select(meta, "restart_V", cid)
      elseif has_variable(meta.footer, pop*"/vg_v")
         # multipop v5 bulk file
         Vbulk = read_variable_select(meta, pop*"/vg_v", cid)
      elseif has_variable(meta.footer, pop*"/V")
         # multipop bulk file
         Vbulk = read_variable_select(meta, pop*"/V", cid)
      else
         # regular bulk file, currently analysator supports pre- and
         # post-multipop files with "V"
         Vbulk = read_variable(meta, "V", cid)
      end

      for f in ("fsaved", "vg_f_saved")
         if has_variable(meta.footer, f) &&
            read_variable_select(meta, f, cid) != 1.0
            @error "VDF not found in the given cell!"
         end
      end
      
      vcellids, vcellf = read_velocity_cells(meta, cid; pop)

      V = get_velocity_cell_coordinates(meta, vcellids; pop)

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
      if has_variable(meta.footer, pop*"/EffectiveSparsityThreshold")
         fThreshold = read_variable_select(meta,
            pop*"/EffectiveSparsityThreshold", cid)
      elseif has_variable(meta.footer, pop*"/vg_effectivesparsitythreshold")
         fThreshold = read_variable_select(meta,
            pop+"/vg_effectivesparsitythreshold", cid)
      else
         verbose && @info "Using a default f threshold value of 1e-16."
         fThreshold = 1e-16
      end

      # Drop all velocity cells which are below the sparsity threshold
      fselect_ = vcellf .≥ fThreshold
      f = vcellf[fselect_]
      V = V[:,fselect_]

      if has_parameter(meta, "t")
         timesim = read_parameter(meta, "t")
         str_title = @sprintf "t= %4.1fs" timesim
      elseif has_parameter(meta, "time")
         timesim = read_parameter(meta, "time")
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
         if has_variable(meta.footer, "B_vol")
            B = read_variable_select(meta, "B_vol", cid)
         elseif has_variable(meta.footer, "vg_b_vol")
            B = read_variable_select(meta, "vg_b_vol", cid)
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

      cb = colorbar(h[4], ax=ax, fraction=0.046, pad=0.04)
      cb_title = cb.ax.set_ylabel("f(v)", fontsize=14)

      if slicetype in ("bperp", "bpar", "bpar1")
         # Draw vector of magnetic field direction
      end
      plt.tight_layout()
   end

end