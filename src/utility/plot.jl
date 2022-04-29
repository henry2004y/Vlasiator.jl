# Plot helpers

"Axis unit type. Currently supported: `SI`, `EARTH`."
@enum AxisUnit SI EARTH
"Color scales type for 2D plots. Currently supported: `Log`, `Linear`, `SymLog`."
@enum ColorScale Log Linear SymLog

"Plotting arguments."
struct PlotArgs
   "Unit of spatial axis"
   axisunit::AxisUnit
   "data array size"
   sizes::Tuple{Int64, Int64}
   "plotting data range"
   plotrange::NTuple{4, Float64}
   "cell IDs in the cut plane"
   idlist::Vector{UInt}
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

"""
    set_args(meta::MetaVLSV, var::String, axisunit::AxisUnit;
       normal::Symbol=:none, origin=0.0)

Set plot-related arguments of `var` in `axisunit`. `normal` and `origin` are used for 2D
slices of 3D data, as specified in [`pcolormeshslice`](@ref).
"""
function set_args(meta::MetaVLSV, var::String, axisunit::AxisUnit; normal::Symbol=:none,
   origin=0.0)

   (;ncells, coordmin, coordmax) = meta

   if normal == :x
      seq = (2,3)
      dir = 1
   elseif normal == :y || (ncells[2] == 1 && ncells[3] != 1) # polar
      seq = (1,3)
      dir = 2
   elseif normal == :z || (ncells[3] == 1 && ncells[2] != 1) # ecliptic
      seq = (1,2)
      dir = 3
   else
      throw(ArgumentError("1D data detected. Please use 1D plot functions."))
   end

   plotrange = (coordmin[seq[1]], coordmax[seq[1]], coordmin[seq[2]], coordmax[seq[2]])
   axislabels = ['X', 'Y', 'Z'][[seq...]]
   # Scale the sizes to the highest refinement level
   sizes = ncells[[seq...]] .<< meta.maxamr # data needs to be refined later

   if normal == :none
      idlist, indexlist = UInt[], Int[]
   else
      idlist, indexlist = let sliceoffset = origin - coordmin[dir]
         getslicecell(meta, sliceoffset, dir, coordmin[dir], coordmax[dir])
      end
   end

   unitstr = axisunit == EARTH ? L"$R_E$" : L"$m$"
   strx = axislabels[1]*" ["*unitstr*"]"
   stry = axislabels[2]*" ["*unitstr*"]"

   str_title = @sprintf "t= %4.1fs" meta.time

   datainfo = readvariablemeta(meta, var)

   cb_title = !isempty(datainfo.variableLaTeX) ?
      datainfo.variableLaTeX * " ["*datainfo.unitLaTeX*"]" : ""

   PlotArgs(axisunit, sizes, plotrange, idlist, indexlist, str_title, strx, stry, cb_title)
end

"""
    set_lim(vmin, vmax, data, colorscale=Linear)

Set colormap limits `vmin`, `vmax` for `data` under scale `colorscale`.
"""
function set_lim(vmin::T, vmax::T, data, colorscale::ColorScale=Linear) where T
   if colorscale in (Linear, SymLog)
      v1 = isinf(vmin) ? minimum(x->isnan(x) ? typemax(T) : convert(T,x), data) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? typemin(T) : convert(T,x), data) : vmax
   else # logarithmic
      datapositive = data[data .> 0.0f0]
      v1 = isinf(vmin) ? minimum(datapositive) : vmin
      v2 = isinf(vmax) ? maximum(x->isnan(x) ? typemin(T) : convert(T,x), data) : vmax
   end

   v1, v2
end

"""
    get_axis(axisunit::AxisUnit, plotrange::NTuple{4,<:Real}, sizes::NTuple{2,<:Integer})
    get_axis(pArgs::PlotArgs)

Return x and y ranges for 2D.
"""
function get_axis(axisunit::AxisUnit, plotrange::NTuple{4,<:Real},
   sizes::NTuple{2,<:Integer})

   if axisunit == EARTH
      x = LinRange(plotrange[1], plotrange[2], sizes[1]) ./ RE
      y = LinRange(plotrange[3], plotrange[4], sizes[2]) ./ RE
   else
      x = LinRange(plotrange[1], plotrange[2], sizes[1])
      y = LinRange(plotrange[3], plotrange[4], sizes[2])
   end
   x, y
end

get_axis(pArgs::PlotArgs) = get_axis(pArgs.axisunit, pArgs.plotrange, pArgs.sizes)

"""
    prep2d(meta::MetaVLSV, var::String, comp=0) -> Array

Obtain data from `meta` of `var` for 2D plotting. Use `comp` to select vector components.
"""
function prep2d(meta::MetaVLSV, var::String, comp=0)
   dataRaw = Vlasiator.getdata2d(meta, var)

   data =
      if ndims(dataRaw) == 3
         if comp in (:x, 1)
            @view dataRaw[1,:,:]
         elseif comp in (:y, 2)
            @view dataRaw[2,:,:]
         elseif comp in (:z, 3)
            @view dataRaw[3,:,:]
         elseif comp in (:mag, 0)
            @views hypot.(dataRaw[1,:,:], dataRaw[2,:,:], dataRaw[3,:,:])
         end
      else
         dataRaw
      end

   data
end

"""
    prep2dslice(meta::MetaVLSV, var::String, normal, comp, pArgs::PlotArgs)

Return `data` of `var` on a uniform 2D mesh on the finest AMR level. Use `normal` to select
the plane orientation, and `comp` to select the component of a vector, same as in
[`pcolormeshslice`](@ref).
"""
function prep2dslice(meta::MetaVLSV, var::String, normal::Symbol, comp, pArgs::PlotArgs)
   (;idlist, indexlist) = pArgs

   data3D = readvariable(meta, var)

   if startswith(var, "fg_") # field quantities, fsgrid
      throw(ArgumentError("FS grid variable $var plotting in cut currently not supported!"))
   else # moments, dccrg grid
      # vlasov grid, AMR
      if ndims(data3D) == 1
         data2D = data3D[indexlist]

         data = refineslice(meta, idlist, data2D, normal)
      elseif ndims(data3D) == 2
         data2D = data3D[:,indexlist]

         if comp in (:x, :y, :z, 1, 2, 3)
            if comp in (:x, 1)
               slice = @view data2D[1,:]
            elseif comp in (:y, 2)
               slice = @view data2D[2,:]
            elseif comp in (:z, 3)
               slice = @view data2D[3,:]
            end
            data = refineslice(meta, idlist, slice, normal)
         elseif comp in (:mag, 0)
            datax = @views refineslice(meta, idlist, data2D[1,:], normal)
            datay = @views refineslice(meta, idlist, data2D[2,:], normal)
            dataz = @views refineslice(meta, idlist, data2D[3,:], normal)
            data = hypot.(datax, datay, dataz)
         end
      end
   end

   data
end

"""
    prep_vdf(meta::MetaVLSV, location::AbstractVector; kwargs...)

Return the cell velocities `v1, v2`, bin ranges `r1, r2`, cell VDF values `fweight`,
and strings of labels and titles for VDF plots.

# Optional arguments
- `unit::AxisUnit`: location unit in `SI`, `EARTH`.
- `unitv::String`: velocity unit in ("km/s", "m/s").
- `limits::Vector{Real}`: velocity space range given in [xmin, xmax, ymin, ymax].
- `slicetype`: symbol for choosing the slice type from `:xy`, `:xz`, `:yz`,
`:bperp`, `:bpar1`, `:bpar2`.
- `center`: symbol for setting the reference frame from `:bulk`, `:peak`.
- `vslicethick`: setting the velocity space slice thickness in the normal direction. If set
to 0, the whole distribution along the normal direction is projected onto a plane. Currently
this is only meaningful when `center` is set such that a range near the bulk/peak normal
velocity is selected!
- `weight::Symbol`: choosing distribution weights from phase space density or particle flux
between `:particle` and `:flux`.
- `flimit`: minimum VDF threshold for plotting.
- `verbose`: display the selection process.
"""
function prep_vdf(meta::MetaVLSV, location::AbstractVector;
   species="proton", unit::AxisUnit=SI, unitv="km/s", slicetype=:default, vslicethick=0.0,
   center=:nothing, weight=:particle, flimit=-1.0, verbose=false)

   ncells = meta.ncells

   if haskey(meta.meshes, species)
      vmesh = meta.meshes[species]
   else
      throw(ArgumentError("Unable to detect population $species"))
   end

   if slicetype ∉ (:default, :xy, :xz, :yz, :bperp, :bpar1, :bpar2)
      throw(ArgumentError("Unknown type $slicetype"))
   end

   unit == EARTH && (location .*= RE)

   # Set unit conversion factor
   unitvfactor = unitv == "km/s" ? 1f3 : 1.0f0

   # Get closest cell ID from input coordinates
   cidReq = getcell(meta, location)
   cidNearest = getnearestcellwithvdf(meta, cidReq)

   # Set normal direction
   if slicetype == :default
      slicetype =
         if ncells[2] == 1 && ncells[3] == 1 # 1D, select xz
            :xz
         elseif ncells[2] == 1 # polar
            :xz
         elseif ncells[3] == 1 # ecliptic
            :xy
         end
   end

   if slicetype in (:xy, :yz, :xz)
      if slicetype == :xy
         dir1, dir2, dir3 = 1, 2, 3
         ŝ = @SVector [0., 0., 1.]
      elseif slicetype == :xz
         dir1, dir2, dir3 = 1, 3, 2
         ŝ = @SVector [0., 1., 0.]
      elseif slicetype == :yz
         dir1, dir2, dir3 = 2, 3, 1
         ŝ = @SVector [1., 0., 0.]
      end
      v1size = vmesh.vblocks[dir1] * vmesh.vblock_size[dir1]
      v2size = vmesh.vblocks[dir2] * vmesh.vblock_size[dir2]
   
      v1min, v1max = vmesh.vmin[dir1], vmesh.vmax[dir1]
      v2min, v2max = vmesh.vmin[dir2], vmesh.vmax[dir2]
   elseif slicetype in (:bperp, :bpar1, :bpar2)
      if hasvariable(meta, "vg_b_vol")
         b̂ = readvariable(meta, "vg_b_vol", cidNearest) |> vec |> normalize!
      elseif hasvariable(meta, "B_vol")
         b̂ = readvariable(meta, "B_vol", cidNearest) |> vec |> normalize!
      end
      if hasvariable(meta, "proton/vg_v")
         v̂ = readvariable(meta, "proton/vg_v", cidNearest) |> vec |> normalize!
      end

      eperp1 = b̂ × v̂
      êperp1 = normalize(eperp1)
      # b̂ ∥ v̂
      norm(eperp1) < 1f-1 && @warn "B is almost parallel to V!"

      êperp2 = normalize(b̂ × êperp1)

      if slicetype == :bperp # slice in b_perp1/b_perp2
         ŝ = SVector{3}(b̂)
         es = SMatrix{3,3}(hcat(êperp1, êperp2, b̂))
      elseif slicetype == :bpar2 # slice in b_parallel/b_perp1 plane
         ŝ = SVector{3}(êperp2)
         es = SMatrix{3,3}(hcat(b̂, êperp1, êperp2))
      else # slice in b_parallel/b_perp2 plane
         ŝ = SVector{3}(êperp1)
         es = SMatrix{3,3}(hcat(êperp2, b̂, êperp1))
      end
      Rinv = let e = @SMatrix Float32[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
         getRotationMatrix(e, es)'
      end
      # heuristic guess
      v1size = vmesh.vblocks[1] * vmesh.vblock_size[1]
      v2size = v1size
      # Compute the ranges along the two perpendicular axis in the slice plane
      r = zeros(Float32, 3, 8)
      r[:,1] = Rinv * Float32[vmesh.vmin[1], vmesh.vmin[2], vmesh.vmin[3]]
      r[:,2] = Rinv * Float32[vmesh.vmin[1], vmesh.vmin[2], vmesh.vmax[3]]
      r[:,3] = Rinv * Float32[vmesh.vmin[1], vmesh.vmax[2], vmesh.vmin[3]]
      r[:,4] = Rinv * Float32[vmesh.vmin[1], vmesh.vmax[2], vmesh.vmax[3]]
      r[:,5] = Rinv * Float32[vmesh.vmax[1], vmesh.vmin[2], vmesh.vmin[3]]
      r[:,6] = Rinv * Float32[vmesh.vmax[1], vmesh.vmin[2], vmesh.vmax[3]]
      r[:,7] = Rinv * Float32[vmesh.vmax[1], vmesh.vmax[2], vmesh.vmin[3]]
      r[:,8] = Rinv * Float32[vmesh.vmax[1], vmesh.vmax[2], vmesh.vmax[3]]

      v1min, v1max = @views extrema(r[1,:])
      v2min, v2max = @views extrema(r[2,:])
   end

   if !isapprox((v1max - v1min) / v1size, (v2max - v2min) / v2size)
      @warn "Noncubic vgrid applied!"
   end
   cellsize = (v1max - v1min) / v1size

   r1 = LinRange(v1min, v1max, v1size+1) / unitvfactor
   r2 = LinRange(v2min, v2max, v2size+1) / unitvfactor

   vcellids, vcellf = readvcells(meta, cidNearest; species)

   V = getvcellcoordinates(meta, vcellids; species)

   if center == :bulk # centered with bulk velocity
      Vbulk =
         if hasvariable(meta, "moments") # From a restart file
            readvariable(meta, "restart_V", cidNearest)
         elseif hasvariable(meta, species*"/vg_v") # Vlasiator 5
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
   elseif slicetype ∈ (:bperp, :bpar1, :bpar2)
      v1select = zeros(eltype(Vselect[1]), length(Vselect))
      v2select = similar(v1select)
      vnormal = similar(v1select)

      for iv in eachindex(Vselect)
         v1select[iv], v2select[iv], vnormal[iv] = Rinv * Vselect[iv]
      end

      if slicetype == :bperp
         strx = L"$v_{B \times V}$"
         stry = L"$v_{B \times (B \times V)}$"
      elseif slicetype == :bpar2
         strx = L"$v_{B}$"
         stry = L"$v_{B \times V}$"
      elseif slicetype == :bpar1
         strx = L"$v_{B \times (B \times V)}$"
         stry = L"$v_{B}$"
      end
   end

   let str = " [$unitv]"
      strx *= str
      stry *= str
   end

   if vslicethick < 0 # Set a proper thickness
      vslicethick =
         if any(ŝ .== 1.0) # Assure that the slice cut through at least 1 vcell
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

   str_title = @sprintf "t = %4.1fs" meta.time

   if verbose
      @info "Original coordinates : $location"
      @info "Original cell        : $(getcellcoordinates(meta, cidReq))"
      @info "Nearest cell with VDF: $(getcellcoordinates(meta, cidNearest))"
      @info "CellID: $cidNearest"

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