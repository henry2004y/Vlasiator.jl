using Vlasiator, UnPack, Printf, StaticArrays
using GLMakie

function Makie.convert_arguments(P::PointBased, meta::MetaVLSV, var)
   data = readvariable(meta, var)
   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   ([Point2f0(i, j) for (i, j) in zip(x, data)],)
end

function Makie.convert_arguments(P::SurfaceLike, meta::MetaVLSV, var;
   axisunit=RE, op=:mag)
   pArgs = Vlasiator.set_args(meta, var, axisunit)
   x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
   data = Vlasiator.prep2d(meta, var, op)
   
   (x, y, data)
end

"""
    viz(meta, var)

Visualize Vlasiator output `var` in `meta` with various options:
* `axisunit`   - unit of axis of type `AxisUnit`
* `colorscale` - scale of colormap of type `ColorScale`
* `normal`     - slice normal direction
* `vmin`       - minimum color value
* `vmax`       - maximum color value
* `op`         - selection of vector components
"""
@recipe(Viz, meta, var) do scene
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
      normal        = :y, # only works in 3D
      vmin          = -Inf,
      vmax          = Inf,
      op            = :mag,
      origin        = 0.0,
   )
end

function Makie.plot!(vlplot::Viz)
   meta        = vlplot[:meta][]
   var         = vlplot[:var][]
   op          = vlplot.op[]
   normal      = vlplot.normal[]
   axisunit    = vlplot.axisunit[]
   vmin        = vlplot.vmin[]
   vmax        = vlplot.vmax[]
   colorscale  = vlplot.colorscale[]
   origin      = vlplot.origin[]

   if ndims(meta) == 1
      data = readvariable(meta, var)

      x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

      lines!(vlplot, x, data)
   elseif ndims(meta) == 2
      pArgs = Vlasiator.set_args(meta, var, axisunit)

      x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
      data = Vlasiator.prep2d(meta, var, op)
   
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

      heatmap!(vlplot, x, y, data, colormap=vlplot.colormap)
   else # 3D
      pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
      x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
      data = Vlasiator.prep2dslice(meta, var, normal, op, pArgs)

      heatmap!(vlplot, x, y, data, colormap=vlplot.colormap)
   end
   vlplot
end

function vlheatmap(meta, var; addcolorbar=true, axisunit=RE, kwargs...)
   pArgs = Vlasiator.set_args(meta, var, axisunit)

   fig = Figure()
   c = viz(fig[1,1], meta, var; axisunit, kwargs...)
   c.axis.title = @sprintf "t= %4.1fs" meta.time
   # TODO: current limitation in Makie 0.15 that no conversion from initial String type
   c.axis.xlabel = pArgs.strx
   c.axis.ylabel = pArgs.stry
   c.axis.autolimitaspect = 1
   if addcolorbar
      cbar = Colorbar(fig[1,2], c.plot, label=pArgs.cb_title_use, tickalign=1)
      colgap!(fig.layout, 7)
   end
   fig
end

"Interactive 2D slice of 3D `var`."
function vlslice(meta, var; normal=:y, axisunit=RE, op=:mag)
   if normal == :x
      dir = 1
   elseif normal == :y
      dir = 2
   else
      dir = 3
   end

   dx = meta.dcoord[dir] / 2^meta.maxamr

   pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin=0.0)
   x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)

   sliceindex = Node(1)
   
   slice = @lift(
      begin
         origin = ($sliceindex-1)*dx + meta.coordmin[dir]
         pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
         data = Vlasiator.prep2dslice(meta, var, normal, op, pArgs)
      end)

   fig = Figure()
   
   heatmap(fig[1, 1], slice)

   nsize = meta.ncells[dir]
   depth = nsize*2^meta.maxamr

   ls = labelslider!(fig, "location:", 2:depth;
      format = x -> "$(x) cells",
      sliderkw = Dict(:startvalue=>depth÷2+1))
   fig[2, 1] = ls.layout

   connect!(sliceindex, ls[1].value)
   
   fig
end