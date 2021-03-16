# Using user recipes from Plots.

using RecipesBase

# Build a recipe which acts on a custom type.
# The function name here is meaningless: it is only used to process a unique 
# set of types early in the pipeline.
# It will work on all functions given the correct dimensions, e.g.
# contourf(data, "rho", xlabel="x")
@recipe function f(meta::MetaData, var::AbstractString; op="mag", axisunit="Re")

   # Check if ecliptic or polar run
   if meta.ycells == 1 && meta.zcells != 1
      plotrange = [meta.xmin, meta.xmax, meta.zmin, meta.zmax]
      sizes = [meta.xcells, meta.ycells]
      PLANE = "XZ"
      axislabels = ['X', 'Z']
   elseif meta.zcells == 1 && meta.ycells != 1
      plotrange = [meta.xmin, meta.xmax, meta.ymin, meta.ymax]
      sizes = [meta.xcells, meta.ycells]
      PLANE = "XY"
      axislabels = ['X', 'Y']
   end

   if var in keys(variables_predefined)
      dataRaw = variables_predefined[var](meta)
   else
      dataRaw = read_variable(meta, var)
   end

   if startswith(var, "fg_") # fsgrid
      data = dataRaw
   else # vlasov grid
      if ndims(dataRaw) == 1 || (ndims(dataRaw) == 2 && size(dataRaw)[1] == 1)       
         data = reshape(dataRaw, sizes[1], sizes[2])
      elseif ndims(dataRaw) == 2
         dataRaw = reshape(dataRaw, 3, sizes...)
         if op == "x"
            data = dataRaw[1,:,:]
         elseif op == "y"
            data = dataRaw[2,:,:]
         elseif op == "z"
            data = dataRaw[3,:,:]
         elseif startswith("mag",op)
            data = hypot.(dataRaw[1,:,:,], dataRaw[2,:,:], dataRaw[3,:,:])
         end
      end
   end

   if axisunit == "Re"
      x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Re
      y = range(plotrange[3], plotrange[4], length=sizes[2]) ./ Re      
   else
      x = range(plotrange[1], plotrange[2], length=sizes[1])
      y = range(plotrange[3], plotrange[4], length=sizes[2])
   end

   @series begin
      seriestype --> :heatmap  # use := if you want to force it
      x, y, data'
   end
end