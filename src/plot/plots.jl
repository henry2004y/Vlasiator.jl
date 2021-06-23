# Using user recipes from Plots.

using RecipesBase

# Build a recipe which acts on a custom type.
@recipe function f(meta::MetaData, var::AbstractString; op=:mag, axisunit=RE)

   # Check if ecliptic or polar run
   if meta.ncells[2] == 1 && meta.ncells[3] != 1
      plotrange = [meta.coordmin[1], meta.coordmax[1], meta.coordmin[3], meta.coordmax[3]]
      sizes = [meta.ncells[1], meta.ncells[2]]
      PLANE = "XZ"
      axislabels = ['X', 'Z']
   elseif meta.ncells[3] == 1 && meta.ncells[2] != 1
      plotrange = [meta.coordmin[1], meta.coordmax[1], meta.coordmin[2], meta.coordmax[2]]
      sizes = [meta.ncells[1], meta.ncells[2]]
      PLANE = "XY"
      axislabels = ['X', 'Y']
   end

   if hasvariable(meta.footer, var)
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
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ Re
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ Re      
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end

   @series begin
      seriestype --> :heatmap  # use := if you want to force it
      x, y, data'
   end
end