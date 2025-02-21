# Type conversion from Vlasiator to Makie

"Conversion for 1D plots"
function Makie.convert_arguments(P::Makie.PointBased, meta::MetaVLSV, var::String)
   data = readvariable(meta, var)
   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   ([Makie.Point2f(i, j) for (i, j) in zip(x, data)],)
end

"Conversion for 2D plots."
function Makie.convert_arguments(P::Makie.CellGrid, meta::MetaVLSV, var::String,
   axisunit::AxisUnit=EARTH, comp::Union{Symbol, Int}=0, normal::Symbol=:y,
   origin::AbstractFloat=0.0)

   if meta.maxamr > 0
      pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
      data = Vlasiator.prep2dslice(meta, var, normal, comp, pArgs)
   else
      pArgs = Vlasiator.set_args(meta, var, axisunit)
      data = Vlasiator.prep2d(meta, var, comp)
   end

   x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)

   (x, y, data)
end

"Conversion for 3D plots."
function Makie.convert_arguments(P::Makie.VolumeLike, meta::MetaVLSV, var::String,
   axisunit::AxisUnit=EARTH, comp::Union{Symbol, Int}=1)
   (; coordmin, coordmax) = meta

   if axisunit == EARTH
      x = (coordmin[1]/Vlasiator.RE, coordmax[1]/Vlasiator.RE)
      y = (coordmin[2]/Vlasiator.RE, coordmax[2]/Vlasiator.RE)
      z = (coordmin[3]/Vlasiator.RE, coordmax[3]/Vlasiator.RE)
   else
      x = (coordmin[1], coordmax[1])
      y = (coordmin[2], coordmax[2])
      z = (coordmin[3], coordmax[3])
   end

   if startswith(var, "fg")
      data = meta[var]
      if comp == 0
         data = [√(data[1,i,j,k]^2 + data[2,i,j,k]^2 + data[3,i,j,k]^2)
            for k in axes(data, 4), j in axes(data, 3), i in axes(data, 2)]
      else
         data = data[comp,:,:,:]
      end
   else
      data, _ = Vlasiator.fillmesh(meta, var)
      # Select on the finest refinement level
      if comp == 0
         data = data[1][end][1,:,:,:]
      else
         data = data[1][end][comp,:,:,:]
      end
   end

   (x, y, z, data)
end