# Plot helpers

"Axis unit type. Currently supported: `SI`, `RE`."
@enum AxisUnit SI RE
"Color scales type for 2D plots. Currently supported: `Log`, `Linear`."
@enum ColorScale Log Linear SymLog

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
function plot_prep2d(meta::MetaVLSV, var, op=:none)
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

   data'
end