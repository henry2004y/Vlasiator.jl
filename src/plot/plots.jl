# Using user recipes from Plots.

using RecipesBase, Printf

# Build a recipe which acts on a custom type.
@recipe function f(meta::MetaVLSV, var::AbstractString; op=:mag, axisunit=RE, normal=:y)
   (;ncells, coordmin, coordmax) = meta
   if ndims(meta) == 1
      data = readvariable(meta, var)

      x = LinRange(coordmin[1], coordmax[1], ncells[1])

      @series begin
         seriestype --> :line  # use := if you want to force it
         x, data
      end
   elseif ndims(meta) == 2
      # Check if ecliptic or polar run
      if ncells[2] == 1 && ncells[3] != 1
         plotrange = (coordmin[1], coordmax[1], coordmin[3], coordmax[3])
         sizes = (ncells[1], ncells[2])
         PLANE = "XZ"
         axislabels = ('X', 'Z')
      elseif ncells[3] == 1 && ncells[2] != 1
         plotrange = (coordmin[1], coordmax[1], coordmin[2], coordmax[2])
         sizes = (ncells[1], ncells[2])
         PLANE = "XY"
         axislabels = ('X', 'Y')
      end

      x, y = Vlasiator.get_axis(axisunit, plotrange, sizes)
      data = Vlasiator.prep2d(meta, var, op)'
      unitstr = axisunit == RE ? "R_E" : "m"

      strx = L"\textrm{%$(axislabels[1])}[%$unitstr]"
      stry = L"\textrm{%$(axislabels[2])}[%$unitstr]"

      @series begin
         seriestype --> :heatmap  # use := if you want to force it
         seriescolor --> :turbo
         xguide --> strx
         yguide --> stry
         title --> @sprintf "t= %4.1fs" meta.time
         x, y, data
      end
   elseif ndims(meta) == 3
      pArgs = set_args(meta, var, axisunit; normal, origin)
      data = prep2dslice(meta, var, normal, pArgs)'
      x, y = Vlasiator.get_axis(axisunit, plotrange, sizes)

      @series begin
         seriestype --> :heatmap  # use := if you want to force it
         seriescolor --> :turbo
         title --> @sprintf "t= %4.1fs" meta.time
         x, y, data
      end
   end
end

@userplot VDFSlice

@recipe function f(h::VDFSlice; species="proton", unit=SI, unitv="km/s", slicetype=:nothing,
   vslicethick=0.0, center=:nothing, fmin=-Inf, fmax=Inf, weight=:particle, flimit=-1.0,
   verbose=false)

   if length(h.args) != 2 || !(typeof(h.args[1]) == MetaVLSV) ||
      !(length(h.args[2]) == 3)
      error("Input argument error for VDF slices.  Got: $(typeof.(h.args))")
   end
   meta, location = h.args

   v1, v2, r1, r2, fweight, strx, stry, str_title =
      Vlasiator.prep_vdf(meta, location;
         species, unit, unitv, slicetype, vslicethick, center, weight, flimit, verbose)

   isinf(fmin) && (fmin = minimum(fweight))
   isinf(fmax) && (fmax = maximum(fweight))

   verbose && @info "Active f range is $fmin, $fmax"

   # setups
   legend := false
   grid := true

   #TODO: Plots 1.22.6 needs improvements on colorbar!
   @series begin
      seriestype := :histogram2d
      seriescolor --> :turbo
      colorbar := true
      colorbar_scale --> :log10
      clims --> (fmin, fmax)
      colorbar_title := "f(v)"
      weights --> fweight
      bins --> r1, r2
      v1, v2
   end

end