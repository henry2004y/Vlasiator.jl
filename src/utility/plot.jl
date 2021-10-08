# Plot helpers

export SI, RE, Log, Linear, SymLog

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