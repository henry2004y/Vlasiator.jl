using Vlasiator, UnPack, Printf
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
   meta = vlplot[:meta][]
   var  = vlplot[:var][]
   
   @unpack ncells, coordmin, coordmax = meta

   if vlplot.normal[] == :x
      sizes = [ncells[2], ncells[3]]
      plotrange = [coordmin[2], coordmax[2], coordmin[3], coordmax[3]]
      sliceoffset = abs(coordmin[1]) + origin

      idlist, indexlist = getslicecell(meta, sliceoffset; xmin=coordmin[1],xmax=coordmax[1])
   elseif vlplot.normal[] == :y
      sizes = [ncells[1], ncells[3]]
      plotrange = [coordmin[1], coordmax[1], coordmin[3], coordmax[3]]
      sliceoffset = abs(coordmin[2]) + origin

      idlist, indexlist = getslicecell(meta, sliceoffset; ymin=coordmin[2],ymax=coordmax[2])
   elseif vlplot.normal[] == :z
      sizes = [ncells[1], ncells[2]]
      plotrange = [coordmin[1], coordmax[1], coordmin[2], coordmax[2]]
      sliceoffset = abs(coordmin[3]) + origin

      idlist, indexlist = getslicecell(meta, sliceoffset; zmin=coordmin[3],zmax=coordmax[3])
   else
      idlist, indexlist = Int64[], Int64[]

      if ncells[2] == 1 && ncells[3] != 1 # polar
         plotrange = [coordmin[1], coordmax[1], coordmin[3], coordmax[3]]
         sizes = [ncells[1], ncells[3]]
      elseif ncells[3] == 1 && ncells[2] != 1 # ecliptic
         plotrange = [coordmin[1], coordmax[1], coordmin[2], coordmax[2]]
         sizes = [ncells[1], ncells[2]]
      else # 1D
         throw(ArgumentError("1D data detected. Please use 1D plot functions."))
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^meta.maxamr # data needs to be refined later

   dataRaw = Vlasiator.getdata2d(meta, var)

   if ndims(dataRaw) == 3
      if vlplot.op[] in (:x, :1)
         data = @view dataRaw[1,:,:]
      elseif vlplot.op[] in (:y, :2)
         data = @view dataRaw[2,:,:]
      elseif vlplot.op[] in (:z, :3)
         data = @view dataRaw[3,:,:]
      elseif vlplot.op[] == :mag
         data = @views hypot.(dataRaw[1,:,:], dataRaw[2,:,:], dataRaw[3,:,:])
      end
   else
      data = dataRaw
   end

   x, y = Vlasiator.get_axis(vlplot.axisunit[], plotrange, sizes)

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

   if vlplot.colorscale[] == Log # Logarithmic plot
      if any(<(0), data)
         throw(DomainError(data, "Nonpositive data detected: use linear scale instead!"))
      end
      datapositive = data[data .> 0.0]
      v1 = isinf(vlplot.vmin[]) ? minimum(datapositive) : vmin
      v2 = isinf(vlplot.vmax[]) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   elseif vlplot.colorscale[] == Linear
      v1 = isinf(vlplot.vmin[]) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vlplot.vmax[]) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
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