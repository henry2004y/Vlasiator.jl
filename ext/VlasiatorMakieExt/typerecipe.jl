# Type conversion from Vlasiator to Makie

"Conversion for 1D plots"
function Makie.convert_arguments(P::Makie.PointBased, meta::MetaVLSV, var::String)
   data = readvariable(meta, var)
   x = LinRange(meta.coordmin[1], meta.coordmax[1], meta.ncells[1])

   ([Makie.Point2f(i, j) for (i, j) in zip(x, data)],)
end

"Conversion for 2D plots."
function Makie.convert_arguments(P::Makie.GridBased, meta::MetaVLSV, var::String,
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
   (;ncells, coordmin, coordmax) = meta

   # Scale the sizes to the highest refinement level
   sizes = ncells .<< meta.maxamr # data needs to be refined later

   if axisunit == EARTH
      x = LinRange(coordmin[1], coordmax[1], sizes[1]) ./ Vlasiator.RE
      y = LinRange(coordmin[2], coordmax[2], sizes[2]) ./ Vlasiator.RE
      z = LinRange(coordmin[3], coordmax[3], sizes[3]) ./ Vlasiator.RE
   else
      x = LinRange(coordmin[1], coordmax[1], sizes[1])
      y = LinRange(coordmin[2], coordmax[2], sizes[2])
      z = LinRange(coordmin[3], coordmax[3], sizes[3])
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