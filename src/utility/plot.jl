# Plot helpers

"Axis unit type. Currently supported: `SI`, `RE`."
@enum AxisUnit SI RE
"Color scales type for 2D plots. Currently supported: `Log`, `Linear`, `SymLog`."
@enum ColorScale Log Linear SymLog

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
   "title"
   str_title::String
   "xlabel"
   strx::String
   "ylabel"
   stry::String
   "colorbar title"
   cb_title::String
end

"Set plot-related arguments of `var` in `axisunit`."
function set_args(meta::MetaVLSV, var, axisunit::AxisUnit; normal::Symbol=:none, origin=0.0)
   (;ncells, coordmin, coordmax) = meta

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
      idlist, indexlist = let sliceoffset = origin - coordmin[dir]
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

   cb_title = !isempty(datainfo.variableLaTeX) ?
      datainfo.variableLaTeX * " ["*datainfo.unitLaTeX*"]" : ""

   PlotArgs(sizes, plotrange, idlist, indexlist, str_title, strx, stry, cb_title)
end

"Set colormap limits for `data`."
function set_lim(vmin, vmax, data, colorscale::ColorScale=Linear)
   if colorscale in (Linear, SymLog)
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   else # logarithmic
      datapositive = data[data .> 0.0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   end

   v1, v2
end

"Return x and y ranges for 2D."
function get_axis(axisunit::AxisUnit, plotrange, sizes)
   if axisunit == RE
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end
   x, y
end

"Obtain data from `meta` of `var` for 2D plotting. Use `op` to select vector components."
function prep2d(meta::MetaVLSV, var, op=:none)
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

   data
end

function prep2dslice(meta::MetaVLSV, var, normal, op, pArgs::PlotArgs)
   (;idlist, indexlist) = pArgs

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
   data
end