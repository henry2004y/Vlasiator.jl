# Using user recipes from Plots.

using RecipesBase, Printf, UnPack

# Build a recipe which acts on a custom type.
@recipe function f(meta::MetaVLSV, var::AbstractString; op=:mag, axisunit=RE)
   @unpack ncells, coordmin, coordmax = meta
   if ndims(meta) == 1
      if hasvariable(meta, var)
         data = readvariable(meta, var)
      else
         data = Vlasiator.variables_predefined[Symbol(var)](meta)
      end

      x = LinRange(coordmin[1], coordmax[1], ncells[1])

      @series begin
         seriestype --> :line  # use := if you want to force it
         x, data
      end
   elseif ndims(meta) == 2
      # Check if ecliptic or polar run
      if ncells[2] == 1 && ncells[3] != 1
         plotrange = [coordmin[1], coordmax[1], coordmin[3], coordmax[3]]
         sizes = [ncells[1], ncells[2]]
         PLANE = "XZ"
         axislabels = ['X', 'Z']
      elseif ncells[3] == 1 && ncells[2] != 1
         plotrange = [coordmin[1], coordmax[1], coordmin[2], coordmax[2]]
         sizes = [ncells[1], ncells[2]]
         PLANE = "XY"
         axislabels = ['X', 'Y']
      end

      if hasvariable(meta, var)
         dataRaw = readvariable(meta, var)
      else
         dataRaw = Vlasiator.variables_predefined[var](meta)
      end

      if startswith(var, "fg_") # fsgrid
         data = dataRaw
      else # vlasov grid
         if ndims(dataRaw) == 1 || (ndims(dataRaw) == 2 && size(dataRaw)[1] == 1)
            data = reshape(dataRaw, sizes[1], sizes[2])
         elseif ndims(dataRaw) == 2
            dataRaw = reshape(dataRaw, 3, sizes...)
            if op == :x
               data = @view dataRaw[1,:,:]
            elseif op == :y
               data = @view dataRaw[2,:,:]
            elseif op == :z
               data = @view dataRaw[3,:,:]
            elseif op == :mag
               data = @views hypot.(dataRaw[1,:,:,], dataRaw[2,:,:], dataRaw[3,:,:])
            end
         end
      end

      if axisunit == RE
         x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Vlasiator.Re
         y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Vlasiator.Re
         unitstr = "R_E"
      else
         x = LinRange(plotrange[1], plotrange[2], sizes[1])
         y = LinRange(plotrange[3], plotrange[4], sizes[2])
         unitstr = "m"
      end

      strx = L"\textrm{%$(axislabels[1])}[%$unitstr]"
      stry = L"\textrm{%$(axislabels[2])}[%$unitstr]"

      @series begin
         seriestype --> :heatmap  # use := if you want to force it
         seriescolor --> :turbo
         xguide --> strx
         yguide --> stry
         title --> @sprintf "t= %4.1fs" meta.time
         x, y, data'
      end
   end
end