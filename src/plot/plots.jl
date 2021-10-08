# Using user recipes from Plots.

using RecipesBase, Printf, UnPack

# Build a recipe which acts on a custom type.
@recipe function f(meta::MetaVLSV, var::AbstractString; op=:mag, axisunit=RE)
   @unpack ncells, coordmin, coordmax = meta
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

      x, y = Vlasiator.get_axis(axisunit, plotrange, sizes)
      unitstr = axisunit == RE ? "R_E" : "m"

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