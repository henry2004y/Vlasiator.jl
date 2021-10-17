using Vlasiator, UnPack, Printf, StaticArrays
using GLMakie

@recipe(VlPlot, meta, var) do scene
   Attributes(
      # generic attributes
      colormap      = Makie.theme(scene, :colormap),
      markersize    = Makie.theme(scene, :markersize),

      # Vlasiator.jl attributes
      elementcolor  = :slategray3,
      boundarycolor = :gray30,
      facetcolor    = :gray30,
      vertexcolor   = :black,
      showboundary  = true,
      variable      = nothing,
   )
end

function Makie.plot!(vlplot::VlPlot)
   meta = vlplot[:meta][]
   var  = vlplot[:var][]

   data = readvariable(meta, var)

   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   lines!(vlplot, x, data)
   vlplot
end

@recipe(VlContourf, meta, var) do scene
   Attributes(;
      # generic attributes
      colormap      = :turbo,
      markersize    = Makie.theme(scene, :markersize),
      colorrange    = Makie.automatic,
      levels        = 5,
      linewidth     = 1.0,
      alpha         = 1.0,

      # Vlasiator.jl attributes
      axisunit      = RE,
      colorscale    = Linear,
      normal        = :none,
      vmin          = -Inf,
      vmax          = Inf,
      op            = :mag,
   )
end

function Makie.plot!(vlplot::VlContourf)
   meta     = vlplot[:meta][]
   var      = vlplot[:var][]
   op       = vlplot.op[]
   normal   = vlplot.normal[]
   axisunit = vlplot.axisunit[]
   vmin     = vlplot.vmin[]
   vmax     = vlplot.vmax[]
   colorscale = vlplot.colorscale[]

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

   if vlplot.normal[] == :none
      idlist, indexlist = Int[], Int[]
   else
      idlist, indexlist = let sliceoffset = abs(coordmin[dir]) + origin
         getslicecell(meta, sliceoffset, dir, coordmin[dir], coordmax[dir])
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^meta.maxamr # data needs to be refined later

   x, y = Vlasiator.get_axis(axisunit, plotrange, sizes)
   data = Vlasiator.plot_prep2d(meta, var, op)'

   if var in ("fg_b", "fg_e", "vg_b_vol", "vg_e_vol") || endswith(var, "vg_v")
      rho_ = findfirst(endswith("rho"), meta.variable)
      if !isnothing(rho_)
         rho = readvariable(meta, meta.variable[rho_])
         rho = reshape(rho, pArgs.sizes[1], pArgs.sizes[2])
         mask = findall(==(0.0), rho)

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

   if colorscale == Log # Logarithmic plot
      if any(<(0), data)
         throw(DomainError(data, "Nonpositive data detected: use linear scale instead!"))
      end
      datapositive = data[data .> 0.0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   elseif colorscale == Linear
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
      nticks = 9
      ticks = range(v1, v2, length=nticks)
   end
   vlplot.colorrange = [v1, v2]

   contourf!(vlplot, x, y, data, colormap=vlplot.colormap)

   vlplot
end


function vl_contourf(meta::MetaVLSV, var; addcolorbar::Bool=true, kwargs...)
   f = Figure()
   c = vlcontourf(f[1,1], meta, var; kwargs...)

   if c.plot.attributes.normal[] == :x
      axislabels = ['Y', 'Z']
   elseif c.plot.attributes.normal[] == :y
      axislabels = ['X', 'Z']
   elseif c.plot.attributes.normal[] == :z
      axislabels = ['X', 'Y']
   else
      if meta.ncells[2] == 1 && meta.ncells[3] != 1 # polar
         PLANE = "XZ"
         axislabels = ['X', 'Z']
      elseif meta.ncells[3] == 1 && meta.ncells[2] != 1 # ecliptic
         PLANE = "XY"
         axislabels = ['X', 'Y']
      end
   end

   unitstr = c.plot.attributes.axisunit[] == RE ? "R_E" : "m"

   c.axis.title = @sprintf "t= %4.1fs" meta.time
   # TODO: current limitation in Makie 0.15 that no conversion from initial String type
   c.axis.xlabel = L"\textrm{%$(axislabels[1])}[%$unitstr]"
   c.axis.ylabel = L"\textrm{%$(axislabels[2])}[%$unitstr]"
   c.axis.autolimitaspect = 1

   datainfo = readvariablemeta(meta, var)
   # TODO: wait for \mathrm to be added in LaTeX engine
   #cb_title = L"%$(datainfo.variableLaTeX) [%$(datainfo.unitLaTeX)]"
   cb_title = "$(datainfo.variableLaTeX) [$(datainfo.unitLaTeX)]"

   if addcolorbar
      Colorbar(f[1, end+1], c.plot, label=cb_title)
   end

   c
end


#=
file = "test/data/bulk.2d.vlsv"
var = "proton/vg_rho"

meta = load(file)

vl_contourf(meta, var)

current_figure()
=#