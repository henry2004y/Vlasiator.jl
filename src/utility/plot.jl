# Plot helpers

"Axis unit type. Currently supported: `SI`, `RE`."
@enum AxisUnit SI RE
"Color scales type for 2D plots. Currently supported: `Log`, `Linear`, `SymLog`."
@enum ColorScale Log Linear SymLog

"Plotting arguments."
struct PlotArgs
   "data array size"
   sizes::Vector{Int}
   "plotting data range"
   plotrange::Vector{Float32}
   "cell IDs in the cut plane"
   idlist::Vector{Int}
   "mapping from original cell order to cut plane"
   indexlist::Vector{Int}
   "title"
   str_title::String
   "xlabel"
   strx::String
   "ylabel"
   stry::String
   "colorbar title"
   cb_title::String
end

"Set plot-related arguments of `var` in `axisunit`."
function set_args(meta::MetaVLSV, var, axisunit::AxisUnit; normal::Symbol=:none, origin=0.0)
   @unpack ncells, coordmin, coordmax = meta

   if normal == :x
      seq = @SVector [2,3]
      dir = 1
   elseif normal == :y || (ncells[2] == 1 && ncells[3] != 1) # polar
      seq = @SVector [1,3]
      dir = 2
   elseif normal == :z || (ncells[3] == 1 && ncells[2] != 1) # ecliptic
      seq = @SVector [1,2]
      dir = 3
   else
      throw(ArgumentError("1D data detected. Please use 1D plot functions."))
   end

   sizes = ncells[seq]
   plotrange = [coordmin[seq[1]], coordmax[seq[1]], coordmin[seq[2]], coordmax[seq[2]]]
   axislabels = ['X', 'Y', 'Z'][seq]

   if normal == :none
      idlist, indexlist = Int[], Int[]
   else
      idlist, indexlist = let sliceoffset = origin - coordmin[dir]
         getslicecell(meta, sliceoffset, dir, coordmin[dir], coordmax[dir])
      end
   end

   # Scale the sizes to the highest refinement level
   sizes *= 2^meta.maxamr # data needs to be refined later

   unitstr = axisunit == RE ? L"$R_E$" : L"$m$"
   strx = axislabels[1]*"["*unitstr*"]"
   stry = axislabels[2]*"["*unitstr*"]"

   str_title = @sprintf "t= %4.1fs" meta.time

   datainfo = readvariablemeta(meta, var)

   cb_title = !isempty(datainfo.variableLaTeX) ?
      datainfo.variableLaTeX * " ["*datainfo.unitLaTeX*"]" : ""

   PlotArgs(sizes, plotrange, idlist, indexlist, str_title, strx, stry, cb_title)
end

"Set colormap limits for `data`."
function set_lim(vmin, vmax, data, colorscale::ColorScale=Linear)
   if colorscale in (Linear, SymLog)
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? +Inf : x, data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   else # logarithmic
      datapositive = data[data .> 0.0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? -Inf : x, data) : vmax
   end

   v1, v2
end

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
function prep2d(meta::MetaVLSV, var, op=:none)
   dataRaw = Vlasiator.getdata2d(meta, var)

   data =
      if ndims(dataRaw) == 3
         if op in (:x, :1)
            @view dataRaw[1,:,:]
         elseif op in (:y, :2)
            @view dataRaw[2,:,:]
         elseif op in (:z, :3)
            @view dataRaw[3,:,:]
         elseif op == :mag
            @views hypot.(dataRaw[1,:,:], dataRaw[2,:,:], dataRaw[3,:,:])
         end
      else
         dataRaw
      end
   data
end

function prep2dslice(meta::MetaVLSV, var, normal, op, pArgs::PlotArgs)
   @unpack idlist, indexlist = pArgs

   data = readvariable(meta, var)

   if startswith(var, "fg_") # field quantities, fsgrid
      throw(ArgumentError("FS grid variable $var plotting in cut currently not supported!"))
   else # moments, dccrg grid
      # vlasov grid, AMR
      if ndims(data) == 1
         data = data[indexlist] # find required cells
      elseif ndims(data) == 2
         data = data[:,indexlist] # find required cells
      end

      # Create the plotting grid
      if ndims(data) == 1
         data = refineslice(meta, idlist, data, normal)
      elseif ndims(data) == 2
         if op in (:x, :y, :z, :1, :2, :3)
            if op in (:x, :1)
               slice = @view data[1,:]
            elseif op in (:y, :2)
               slice = @view data[2,:]
            elseif op in (:z, :3)
               slice = @view data[3,:]
            end
            data = refineslice(meta, idlist, slice, normal)
         elseif op == :mag
            datax = @views refineslice(meta, idlist, data[1,:], normal)
            datay = @views refineslice(meta, idlist, data[2,:], normal)
            dataz = @views refineslice(meta, idlist, data[3,:], normal)
            data = hypot.(datax, datay, dataz)
         end

      elseif ndims(data) == 3
         @error "not implemented yet!"
      end
   end
   data
end

"""
    prep_vdf(meta::MetaVLSV, location; kwargs...)

Return the cell velocities `v1, v2`, bin ranges `r1, r2`, cell VDF values `fweight`,
and strings of labels and titles for VDF plots.

# Optional arguments
- `unit::AxisUnit`: location unit in `SI`, `RE`.
- `unitv::String`: velocity unit in ("km/s", "m/s").
- `limits::Vector{Real}`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype`: symbol for choosing the slice type from :xy, :xz, :yz, :bperp, :bpar, :bpar1.
- `center`: symbol for setting the reference frame from :bulk, :peak.
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected!
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
- `flimit`: minimum VDF threshold for plotting.
- `verbose`: display the selection process.
"""
function prep_vdf(meta::MetaVLSV, location;
   species="proton", unit::AxisUnit=SI, unitv="km/s", slicetype=:nothing, vslicethick=0.0,
   center=:nothing, weight=:particle, flimit=-1.0, verbose=false)

   ncells = meta.ncells

   if haskey(meta.meshes, species)
      vmesh = meta.meshes[species]
   else
      throw(ArgumentError("Unable to detect population $species"))
   end

   unit == RE && (location .*= Re)

   # Set unit conversion factor
   unitvfactor = unitv == "km/s" ? 1f3 : 1.0f0

   # Get closest cell ID from input coordinates
   cidReq = getcell(meta, location)
   cidNearest = getnearestcellwithvdf(meta, cidReq)

   # Set normal direction
   if slicetype == :nothing
      slicetype =
         if ncells[2] == 1 && ncells[3] == 1 # 1D, select xz
            :xz
         elseif ncells[2] == 1 # polar
            :xz
         elseif ncells[3] == 1 # ecliptic
            :xy
         end
   end

   if slicetype == :xy
      dir1, dir2, dir3 = 1, 2, 3
      sliceNormal = [0., 0., 1.]
   elseif slicetype == :xz
      dir1, dir2, dir3 = 1, 3, 2
      sliceNormal = [0., 1., 0.]
   elseif slicetype == :yz
      dir1, dir2, dir3 = 2, 3, 1
      sliceNormal = [1., 0., 0.]
   elseif slicetype in (:bperp, :bpar, :bpar1)
      if hasvariable(meta, "B_vol")
         B = readvariable(meta, "B_vol", cidNearest)
      elseif hasvariable(meta, "vg_b_vol")
         B = readvariable(meta, "vg_b_vol", cidNearest)
      end
      BxV = B × Vbulk
      if slicetype == :bperp # slice in b_perp1/b_perp2
         sliceNormal = B ./ norm(B)
      elseif slicetype == :bpar1 # slice in b_parallel/b_perp1 plane
         sliceNormal = B × BxV
         sliceNormal ./= norm(sliceNormal)
      else # slice in b_parallel/b_perp2 plane
         sliceNormal = BxV ./ norm(BxV)
      end
   end

   v1size = vmesh.vblocks[dir1] * vmesh.vblock_size[dir1]
   v2size = vmesh.vblocks[dir2] * vmesh.vblock_size[dir2]

   v1min, v1max = vmesh.vmin[dir1], vmesh.vmax[dir1]
   v2min, v2max = vmesh.vmin[dir2], vmesh.vmax[dir2]

   @assert (v1max - v1min) / v1size ≈ (v2max - v2min) / v2size "Noncubic vgrid detected!"
   cellsize = (v1max - v1min) / v1size

   vcellids, vcellf = readvcells(meta, cidNearest; species)

   V = getvcellcoordinates(meta, vcellids; species)

   if center == :bulk # centered with bulk velocity
      Vbulk =
         if hasvariable(meta, "moments")
            # From a restart file
            readvariable(meta, "restart_V", cidNearest)
         elseif hasvariable(meta, species*"/vg_v")
            # Vlasiator 5
            readvariable(meta, species*"/vg_v", cidNearest)
         elseif hasvariable(meta, species*"/V")
            readvariable(meta, species*"/V", cidNearest)
         else
            readvariable(meta, "V", cidNearest)
         end

      for i in eachindex(V)
         V[i] = V[i] .- Vbulk
      end
   elseif center == :peak # centered on highest VDF-value
      Vpeak = maximum(vcellf)
      for i in eachindex(V)
         V[i] = V[i] .- Vpeak
      end
   end

   # Set sparsity threshold
   if flimit < 0
      flimit =
         if hasvariable(meta, species*"/vg_effectivesparsitythreshold")
            readvariable(meta, species*"/vg_effectivesparsitythreshold", cidNearest)
         elseif hasvariable(meta, species*"/EffectiveSparsityThreshold")
            readvariable(meta, species*"/EffectiveSparsityThreshold", cidNearest)
         else
            1e-16
         end
   end

   # Drop velocity cells which are below the sparsity threshold
   findex_ = vcellf .≥ flimit
   fselect = vcellf[findex_]
   Vselect = V[findex_]

   if slicetype ∈ (:xy, :yz, :xz)
      v1select = [v[dir1] for v in Vselect]
      v2select = [v[dir2] for v in Vselect]
      vnormal = [v[dir3] for v in Vselect]

      strx, stry = getindex(["vx", "vy", "vz"], [dir1, dir2])
   elseif slicetype ∈ (:Bperp, :Bpar, :Bpar1)
      #TODO: NOT working yet!
      if slicetype == :Bperp
         v1select = Vrot2[1,:] # the X axis of the slice is BcrossV=perp1
         v2select = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[3,:] # the Z axis of the slice is B
         strx = L"$v_{B \times V}$"
         stry = L"$v_{B \times (B \times V)}$"
      elseif slicetype == :Bpar
         v1select = Vrot2[3,:] # the X axis of the slice is B
         v2select = Vrot2[2,:] # the Y axis of the slice is Bcross(BcrossV)=perp2
         vnormal = Vrot2[1,:] # the Z axis of the slice is -BcrossV=perp1
         strx = L"$v_{B}$"
         stry = L"$v_{B \times V}$"
      elseif slicetype == :Bpar1
         v1select = Vrot2[3,:] # the X axis of the slice is B
         v2select = Vrot2[1,:] # the Y axis of the slice is BcrossV=perp1
         vnormal = Vrot2[2,:] # the Z axis of the slice is Bcross(BcrossV)=perp2
         strx = L"$v_{B}$"
         stry = L"$v_{B \times (B \times V)}$"
      end
   end

   let str = " [$unitv]"
      strx *= str
      stry *= str
   end

   if vslicethick < 0 # Set a proper thickness
      vslicethick =
         if any(sliceNormal .== 1.0) # Assure that the slice cut through at least 1 vcell
            cellsize
         else # Assume cubic vspace grid, add extra space
            cellsize*(√3+0.05)
      end
   end

   # Weights using particle flux or phase-space density
   fweight = weight == :flux ? fselect*norm([v1select, v2select, vnormal]) : fselect

   # Select cells within the slice area
   if vslicethick > 0.0
      v1, v2, fweight = let ind_ = @. (abs(vnormal) ≤ 0.5*vslicethick)
         v1select[ind_], v2select[ind_], fweight[ind_]
      end
   else
      v1, v2 = v1select, v2select
   end
   
   v1 ./= unitvfactor
   v2 ./= unitvfactor

   r1 = LinRange(v1min, v1max, v1size+1) / unitvfactor
   r2 = LinRange(v2min, v2max, v2size+1) / unitvfactor

   str_title = @sprintf "t = %4.1fs" meta.time

   if verbose
      @info "Original coordinates : $location"
      @info "Original cell        : $(getcellcoordinates(meta, cidReq))"
      @info "Nearest cell with VDF: $(getcellcoordinates(meta, cidNearest))"
      let
         x, y, z = getcellcoordinates(meta, cidNearest)
         @info "cellid $cidNearest, x = $x, y = $y, z = $z"
      end

      if center == :bulk
         @info "Transforming to plasma frame, travelling at speed $Vbulk"
      elseif center == :peak
         @info "Transforming to peak f-value frame, travelling at speed $Vpeak"
      end

      @info "Using VDF threshold value of $flimit."

      if vslicethick > 0
         @info "Performing slice with a counting thickness of $vslicethick"
      else
         @info "Projecting total VDF to a single plane"
      end
   end

   v1, v2, r1, r2, fweight, strx, stry, str_title
end