# Interactive plots with Observables

"""
    vlslice(meta, var; normal=:y, axisunit=SI, comp=0)

Interactive 2D slice of 3D `var` in `normal` direction.
"""
function vlslice(meta::MetaVLSV, var::String;
   normal::Symbol=:y, axisunit::AxisUnit=SI, comp::Union{Symbol, Int}=0)
   dir, str1, str2 =
      if normal == :x
         1, "y", "z"
      elseif normal == :y
         2, "x", "z"
      else
         3, "x", "y"
      end

   unitx = axisunit == EARTH ? " [Re]" : " [m]"

   dx = meta.dcoord[dir] / 2^meta.maxamr

   pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin=0.0)
   x, y = Vlasiator.get_axis(axisunit, pArgs.plotrange, pArgs.sizes)

   nsize = meta.ncells[dir]
   depth = nsize*2^meta.maxamr

   fig = Makie.Figure()
   ax = Makie.Axis(fig[1, 1], aspect=Makie.DataAspect())
   ax.xlabel = str1*unitx
   ax.ylabel = str2*unitx

   lsgrid = Makie.SliderGrid(fig[2, 1],
      (label = "location in normal direction $(String(normal))",
       range=1:depth,
       format=x -> "$(x) cells"),
   )

   sliderobservables = [s.value for s in lsgrid.sliders]

   slice = Makie.lift(sliderobservables...) do slvalues...
      begin
         origin = (slvalues[1]-1)*dx + meta.coordmin[dir]
         pArgs = Vlasiator.set_args(meta, var, axisunit; normal, origin)
         Vlasiator.prep2dslice(meta, var, normal, comp, pArgs)
      end
   end

   Makie.heatmap!(ax, x, y, slice, colormap=:turbo)

   Makie.set_close_to!(lsgrid.sliders[1], .5depth)

   fig, ax
end

"""
    vdfslices(meta, location; fmin=1f-16, species="proton", unit=SI, verbose=false)

Three orthogonal slices of VDFs from `meta` at `location`.
# Optional Arguments
- `fmin`: minimum VDF threshold for plotting.
- `species`: name of particle.
- `unit`: unit of input `location`, `SI` or `EARTH`.
"""
function vdfslices(meta::MetaVLSV, location::AbstractVector; fmin::AbstractFloat=1f-16,
   species::String="proton", unit::AxisUnit=SI, verbose::Bool=false)
   if haskey(meta.meshes, species)
      vmesh = meta.meshes[species]
   else
      throw(ArgumentError("Unable to detect population $species"))
   end

   unit == EARTH && (location .*= Vlasiator.RE)

   # Calculate cell ID from given coordinates
   cidReq = getcell(meta, location)
   cidNearest = getnearestcellwithvdf(meta, cidReq, species)

   cellused = getcellcoordinates(meta, cidNearest)

   if verbose
      @info "Original coordinates : $location"
      @info "Original cell        : $(getcellcoordinates(meta, cidReq))"
      @info "Nearest cell with VDF: $cellused"
      let
         x, y, z = getcellcoordinates(meta, cidNearest)
         @info "cellid $cidNearest, x = $x, y = $y, z = $z"
      end
   end

   vcellids, vcellf = readvcells(meta, cidNearest; species)

   f = Vlasiator.reconstruct(vmesh, vcellids, vcellf)

   fig = Makie.Figure()
   ax = Makie.Axis3(fig[1, 1], aspect=(1,1,1), title = "VDF at $cellused in log scale")
   ax.xlabel = "vx [m/s]"
   ax.ylabel = "vy [m/s]"
   ax.zlabel = "vz [m/s]"

   x = LinRange(vmesh.vmin[1], vmesh.vmax[1], vmesh.vblocksize[1]*vmesh.vblocks[1])
   y = LinRange(vmesh.vmin[2], vmesh.vmax[2], vmesh.vblocksize[2]*vmesh.vblocks[2])
   z = LinRange(vmesh.vmin[3], vmesh.vmax[3], vmesh.vblocksize[3]*vmesh.vblocks[3])

   lsgrid = Makie.SliderGrid(fig[2, 1],
      (label="vx", range=1:length(x), format=i -> "$(round(x[i], digits=2))"),
      (label="vy", range=1:length(y), format=i -> "$(round(y[i], digits=2))"),
      (label="vz", range=1:length(z), format=i -> "$(round(z[i], digits=2))"),
   )

   for i in eachindex(f)
      if f[i] < fmin f[i] = fmin end
   end
   data = log10.(f)

   plt = Makie.volumeslices!(ax, x, y, z, data, colormap=:viridis)
   #TODO: improve on colormap setup!
   cbar = Makie.Colorbar(fig, plt,
      label="f(v)",
      minorticksvisible=true)

   fig[1, 2] = cbar

   # connect sliders to volumeslices update methods
   sl_yz, sl_xz, sl_xy = lsgrid.sliders

   Makie.on(sl_yz.value) do v; plt[:update_yz][](v) end
   Makie.on(sl_xz.value) do v; plt[:update_xz][](v) end
   Makie.on(sl_xy.value) do v; plt[:update_xy][](v) end

   Makie.set_close_to!(sl_yz, .5length(x))
   Makie.set_close_to!(sl_xz, .5length(y))
   Makie.set_close_to!(sl_xy, .5length(z))

   fig
end