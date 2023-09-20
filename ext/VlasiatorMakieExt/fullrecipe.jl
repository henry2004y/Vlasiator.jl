# Full recipes for customized Vlasiator plots.

Makie.@recipe(Viz, meta, var) do scene
   Makie.Attributes(;
      # generic attributes
      colormap      = :turbo,
      markersize    = Makie.theme(scene, :markersize),
      colorrange    = Makie.automatic,
      levels        = 5,
      linewidth     = 1.0,
      alpha         = 1.0,

      # Vlasiator.jl attributes
      axisunit      = EARTH,
      colorscale    = Linear,
      normal        = :y, # only works in 3D
      vmin          = -Inf,
      vmax          = Inf,
      comp          = 0,
      origin        = 0.0,
   )
end

function Makie.plot!(vlplot::Viz)
   meta        = vlplot[:meta][]
   var         = vlplot[:var][]
   comp        = vlplot.comp[]
   normal      = vlplot.normal[]
   axisunit    = vlplot.axisunit[]
   vmin        = vlplot.vmin[]
   vmax        = vlplot.vmax[]
   colorscale  = vlplot.colorscale[]
   origin      = vlplot.origin[]

   if ndims(meta) == 1
      data = readvariable(meta, var)

      x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

      Makie.lines!(vlplot, x, data)
   elseif ndims(meta) == 2
      pArgs = Vlasiator.set_args(meta, var, axisunit)

      x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
      data = Vlasiator.prep2d(meta, var, comp)

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

      Makie.heatmap!(vlplot, x, y, data, colormap=vlplot.colormap)
   else # 3D
      pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
      x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)
      data = Vlasiator.prep2dslice(meta, var, normal, comp, pArgs)

      Makie.heatmap!(vlplot, x, y, data, colormap=vlplot.colormap)
   end
   vlplot
end

function vlslices(meta::MetaVLSV, var::String; fig=nothing, axisunit::AxisUnit=SI,
   comp::Union{Symbol, Int}=0, origin::AbstractVector=[0.0, 0.0, 0.0],
   addcolorbar::Bool=false, colorscale::ColorScale=Linear, vmin::Real=-Inf, vmax::Real=Inf)
   if axisunit == EARTH
      unitx = " [Re]"
      origin .*= Vlasiator.RE
   else
      unitx = " [m]"
   end

   pArgs1 = Vlasiator.set_args(meta, var, axisunit; normal=:x, origin=origin[1])
   pArgs2 = Vlasiator.set_args(meta, var, axisunit; normal=:y, origin=origin[2])
   pArgs3 = Vlasiator.set_args(meta, var, axisunit; normal=:z, origin=origin[3])

   x, y = Vlasiator.get_axis(axisunit, pArgs3.plotrange, pArgs3.sizes)
   x, z = Vlasiator.get_axis(axisunit, pArgs2.plotrange, pArgs2.sizes)

   d1 = Vlasiator.prep2dslice(meta, var, :x, comp, pArgs1)
   d2 = Vlasiator.prep2dslice(meta, var, :y, comp, pArgs2)
   d3 = Vlasiator.prep2dslice(meta, var, :z, comp, pArgs3)

   isnothing(fig) && (fig = Makie.Figure(fontsize=22))
   ax = Makie.Axis3(fig[1,1], aspect=(1, 1, 1), elevation=pi/6, perspectiveness=0.5)

   ax.xlabel = "x"*unitx
   ax.ylabel = "y"*unitx
   ax.zlabel = "z"*unitx

   Makie.xlims!(ax, x[1], x[end])
   Makie.ylims!(ax, y[1], y[end])
   Makie.zlims!(ax, z[1], z[end])

   vmin1, vmax1 = Vlasiator.set_lim(vmin, vmax, d1, colorscale)
   vmin2, vmax2 = Vlasiator.set_lim(vmin, vmax, d2, colorscale)
   vmin3, vmax3 = Vlasiator.set_lim(vmin, vmax, d3, colorscale)
   colormap = :turbo
   colorrange = (min(vmin1, vmin2, vmin3), max(vmax1, vmax2, vmax3))

   h1 = Makie.heatmap!(ax, y, z, d1; colormap, colorrange, transformation=(:yz, origin[1]))
   h2 = Makie.heatmap!(ax, x, z, d2; colormap, colorrange, transformation=(:xz, origin[2]))
   h3 = Makie.heatmap!(ax, x, y, d3; colormap, colorrange, transformation=(:xy, origin[3]))

   if addcolorbar
      cbar = Makie.Colorbar(fig[1,2], h3, label=latexstring(pArgs1.cb_title), tickalign=1)
      Makie.colgap!(fig.layout, 7)
   end

   fig, ax
end

function vdfvolume(meta::MetaVLSV, location::AbstractVector; species::String="proton",
   unit::AxisUnit=SI, flimit::Real=-1.0, verbose::Bool=false, fig=nothing)

   if haskey(meta.meshes, species)
      vmesh = meta.meshes[species]
   else
      throw(ArgumentError("Unable to detect population $species"))
   end

   loc = unit == EARTH ?
      location .* Vlasiator.RE :
      location

   # Calculate cell ID from given coordinates
   cidReq = getcell(meta, loc)
   cidNearest = getnearestcellwithvdf(meta, cidReq, species)
   ccoords = getcellcoordinates(meta, cidNearest)

   if verbose
      @info "Species             : $species"
      @info "Original coordinates: $loc"
      @info "Original cell       : $cidReq"
      @info "Actual cell         : $cidNearest"
      @info "Actual coordinates  : $ccoords"
   end

   vcellids, vcellf = readvcells(meta, cidNearest; species)

   V = getvcellcoordinates(meta, vcellids; species)

   # Set sparsity threshold
   if flimit < 0
      flimit =
         if hasvariable(meta, species*"/vg_effectivesparsitythreshold")
            readvariable(meta, species*"/vg_effectivesparsitythreshold", cidNearest)
         elseif hasvariable(meta, species*"/EffectiveSparsityThreshold")
            readvariable(meta, species*"/EffectiveSparsityThreshold", cidNearest)
         else
            1f-16
         end
   end

   # Drop velocity cells which are below the sparsity threshold
   findex_ = vcellf .≥ flimit
   fselect = vcellf[findex_]
   Vselect = V[findex_]

   cmap = Makie.resample_cmap(:turbo, 101; alpha=(0, 1))

   isnothing(fig) && (fig = Makie.Figure(fontsize=22))
   if unit == SI
      ax = Makie.Axis3(fig[1, 1], aspect=(1,1,1), title="VDF at $ccoords (m) in log scale",
         titlesize=26)
   elseif unit == EARTH
      coords = round.(ccoords ./ Vlasiator.RE, digits=1) 
      ax = Makie.Axis3(fig[1, 1], aspect=(1,1,1), title="VDF at $coords (RE) in log scale",
         titlesize=26)
   end
   ax.xlabel = "vx [m/s]"
   ax.ylabel = "vy [m/s]"
   ax.zlabel = "vz [m/s]"

   plt = Makie.meshscatter!(ax, Vselect, color=log10.(fselect),
      marker=Makie.Rect3f(Makie.Vec3f(0), Makie.Vec3f(4*vmesh.dv[1])),
      colormap=cmap,
      transparency=true, shading=false)

   cbar = Makie.Colorbar(fig, plt, label="f(v)")

   fig[1, 2] = cbar

   fig, ax
end

struct LogMinorTicks end

function Makie.get_minor_tickvalues(::LogMinorTicks, scale, tickvalues, vmin, vmax)
   vals = Float64[]
   extended_tickvalues = [
      tickvalues[1] - (tickvalues[2] - tickvalues[1]);
      tickvalues;
      tickvalues[end] + (tickvalues[end] - tickvalues[end-1]);
   ]
   for (lo, hi) in zip(
        @view(extended_tickvalues[1:end-1]),
        @view(extended_tickvalues[2:end])
     )
       interval = hi - lo
       steps = log10.(LinRange(10^lo, 10^hi, 11))
       append!(vals, steps[2:end-1])
   end
   return filter(x -> vmin < x < vmax, vals)
end

custom_formatter(values) = map(
  v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
  values
)

function vdfslice(meta, location;
   species="proton", unit=SI, unitv="km/s", slicetype=:default, vslicethick=0.0,
   center=:nothing, vmin=-Inf, vmax=Inf, weight=:particle, flimit=-1.0, verbose=false,
   fig=nothing)

   v1, v2, r1, r2, fweight, strx, stry, str_title =
      Vlasiator.prep_vdf(meta, location;
         species, unit, unitv, slicetype, vslicethick, center, weight, flimit, verbose)

   isinf(vmin) && (vmin = minimum(fweight))
   isinf(vmax) && (vmax = maximum(fweight))

   verbose && @info "Active f range is $vmin, $vmax"

   h = fit(Histogram, (v1, v2), weights(fweight), (r1, r2))

   clims = (vmin, maximum(h.weights))

   data = [isinf(x) ? NaN : x for x in log10.(h.weights)]

   isnothing(fig) && (fig = Makie.Figure())

   ax, hm = Makie.heatmap(fig[1, 1], r1, r2, data;
      colormap=:turbo, colorrange=log10.(clims))

   cb = Makie.Colorbar(fig[1, 2], hm;
      label="f(v)",
      tickformat=custom_formatter,
      minorticksvisible=true,
      minorticks=LogMinorTicks() )

   ax.title = str_title
   ax.xlabel = strx
   ax.ylabel = stry

   fig, ax
end