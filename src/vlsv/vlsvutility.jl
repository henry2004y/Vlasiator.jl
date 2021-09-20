# Utility functions for processing VLSV data.

using LinearAlgebra: dot
using WriteVTK, Printf
using LazyGrids: ndgrid

const qₑ = -1.60217662e-19  # electron charge, [C]
const mₑ = 9.10938356e-31   # electron mass, [kg]
const qᵢ = 1.60217662e-19   # proton mass, [C]
const mᵢ = 1.673557546e-27  # proton mass, [kg]
const c  = 3e8              # speed of light, [m/s]
const μ₀ = 4π*1e-7          # Vacuum permeability, [H/m]
const kB = 1.38064852e-23   # Boltzmann constant, [m²kg/(s²K)]
const Re = 6.371e6          # Earth radius, [m]

export getcell, getslicecell, getlevel, refineslice, getcellcoordinates,
   getchildren, getparent, isparent, getsiblings,
   getcellinline, getnearestcellwithvdf, getcellwithvdf, write_vtk, issame

"""
    getcell(meta, location) -> Int

Return cell ID containing the given spatial `location`, excluding domain boundaries.
Only accept 3D location.
"""
function getcell(meta::MetaVLSV, loc)
   @unpack coordmin, coordmax, dcoord, ncells, cellid, maxamr = meta

   @assert coordmin[1] < loc[1] < coordmax[1] "x coordinate out of bound!"
   @assert coordmin[2] < loc[2] < coordmax[2] "y coordinate out of bound!"
   @assert coordmin[3] < loc[3] < coordmax[3] "z coordinate out of bound!"

   dx, dy, dz = dcoord

   # Get cell indices
   indices = @inbounds floor.(Int,
      [(loc[1] - coordmin[1])/dx,
       (loc[2] - coordmin[2])/dy,
       (loc[3] - coordmin[3])/dz])
   # Get the cell id
   cid = @inbounds indices[1] + indices[2]*ncells[1] + indices[3]*ncells[1]*ncells[2] + 1

   # Going through AMR levels as needed
   ilevel = 0
   ncells_lowerlevel = 0
   ncell = prod(ncells)

   @inbounds while ilevel < maxamr + 1
      if cid in cellid
         break
      else
         ncells_lowerlevel += 2^(3*ilevel)*ncell
         ilevel += 1
         dx *= 0.5; dy *= 0.5; dz *= 0.5

         indices = floor.(Int,
            [(loc[1] - coordmin[1])/dx,
             (loc[2] - coordmin[2])/dy,
             (loc[3] - coordmin[3])/dz ])

         cid = ncells_lowerlevel + indices[1] +
            2^(ilevel)*xcells*indices[2] +
            4^(ilevel)*xcells*ycells*indices[3] + 1
      end
   end

   cid
end

"""
    getlevel(meta, cid) -> Int

Return the AMR level of a given cell ID. Note that this function does not check if the VLSV
file of `meta` actually contains `cid`: it may be shadowed by refined children.
"""
function getlevel(meta::MetaVLSV, cid::Integer)
   ncell = prod(meta.ncells)
   ilevel = 0
   while cid > 0
      cid -= 2^(3*ilevel)*ncell
      ilevel += 1
   end
   ilevel - 1
end

"""
    getparent(meta, cid) -> Int

Return the parent cell ID of given child `cid`.
"""
function getparent(meta::MetaVLSV, cid::Integer)
   @inbounds xcell, ycell = meta.ncells[1], meta.ncells[2]
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)
   parentlvl = mylvl - 1

   if parentlvl < 0
      throw(ArgumentError("Cell ID $cid has no parent cell!"))
   else
      # get the first cellid on my level
      cid1st = get1stcell(mylvl, ncell) + 1
      # get row and column sequence on my level (starting with 0)
      xcell = xcell*2^mylvl
      ycell = ycell*2^mylvl

      myseq = cid - cid1st
      ix = myseq % xcell
      iz = myseq ÷ (xcell*ycell)
      iy = (myseq - iz*xcell*ycell) ÷ xcell
      # indexes on the parent level
      ixparent = ix ÷ 2
      iyparent = iy ÷ 2
      izparent = iz ÷ 2

      # get the first cellid on parent level
      cid1st -= ncell*8^parentlvl
      # get parent cellid (may not exist!!!)
      parentid = cid1st + izparent*xcell*ycell÷4 + iyparent*xcell÷2 + ixparent
   end
   parentid
end

"""
    getchildren(meta, cid) -> Vector{Int}

Return direct children of `cid`.
"""
function getchildren(meta::MetaVLSV, cid::Integer)
   xcell, ycell, zcell = meta.ncells
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)

   # get the first cell ID on the my level
   cid1st = 1
   for i = 0:mylvl-1
      cid1st += ncell*8^i
   end
   # get my row and column sequence on my level (starting with 0)
   xcell = xcell*2^mylvl
   ycell = ycell*2^mylvl

   myseq = cid - cid1st
   ix = myseq % xcell
   iz = myseq ÷ (xcell*ycell)
   iy = (myseq - iz*xcell*ycell) ÷ xcell

   # get the children sequences on the finer level
   ix *= 2
   iy *= 2
   iz *= 2

   nchildren = 2^ndims(meta)
   cid = @MVector zeros(Int, nchildren)
   # get the first cell ID on the finer level
   cid1st += ncell*8^mylvl
   ix_, iy_ = [ix, ix+1], [iy, iy+1]
   iz_ = zcell != 1 ? [iz, iz+1] : [iz]
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*xcell*ycell*4 + i[2]*xcell*2 + i[1]
   end
   SVector{nchildren}(cid)
end

"""
    getsiblings(meta, cid) -> Vector{Int}

Return sibling cells of a given `cid`, including itself.
"""
function getsiblings(meta::MetaVLSV, cid::Integer)
   xcell, ycell, zcell = meta.ncells
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)

   mylvl == 0 && throw(ArgumentError("CellID $cid is not a child cell!"))

   xcell = xcell * 2^mylvl
   ycell = ycell * 2^mylvl

   # 1st cellid on my level
   cid1st = get1stcell(mylvl, ncell) + 1

   # xyz sequences on my level (starting with 0)
   myseq = cid - cid1st
   ix = myseq % xcell
   iz = myseq ÷ (xcell*ycell)
   iy = (myseq - iz*xcell*ycell) ÷ xcell

   ix1 = iseven(ix) ? ix + 1 : ix - 1
   iy1 = iseven(iy) ? iy + 1 : iy - 1
   iz1 = iseven(iz) ? iz + 1 : iz - 1
   # reorder
   ix, ix1 = minmax(ix, ix1)
   iy, iy1 = minmax(iy, iy1)
   iz, iz1 = minmax(iz, iz1)

   nsiblings = 2^ndims(meta)
   cid = @MVector zeros(Int, nsiblings)
   ix_, iy_ = [ix, ix1], [iy, iy1]
   iz_ = zcell != 1 ? [iz, iz1] : [iz]
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*xcell*ycell + i[2]*xcell + i[1]
   end
   SVector{nsiblings}(cid)
end

"""
    isparent(meta, cid) -> Bool

Check if `cid` is a parent cell.
"""
function isparent(meta::MetaVLSV, cid::Integer)
   ncell_accum = get1stcell(meta.maxamr, prod(meta.ncells))

   cid ∉ meta.cellid && 0 < cid ≤ ncell_accum
end

"""
    getcellcoordinates(meta, cid) -> Vector{Float}

Return a given cell's coordinates.
"""
function getcellcoordinates(meta::MetaVLSV, cid::Integer)
   @unpack ncells, coordmin, coordmax = meta
   cid -= 1 # for easy divisions

   xcell, ycell, zcell = ncells
   reflevel = 0
   subtraction = prod(ncells) * (2^reflevel)^3
   # sizes on the finest level
   while cid ≥ subtraction
      cid -= subtraction
      reflevel += 1
      subtraction *= 8
      xcell *= 2
      ycell *= 2
      zcell *= 2
   end

   indices = @SVector [
      cid % xcell,
      cid ÷ xcell % ycell,
      cid ÷ (xcell*ycell) ]

   coords = @inbounds @SVector [
      coordmin[1] + (indices[1] + 0.5) * (coordmax[1] - coordmin[1]) / xcell,
      coordmin[2] + (indices[2] + 0.5) * (coordmax[2] - coordmin[2]) / ycell,
      coordmin[3] + (indices[3] + 0.5) * (coordmax[3] - coordmin[3]) / zcell ]

   coords
end

function isInsideDomain(meta::MetaVLSV, point)
   @unpack coordmin, coordmax = meta

   if coordmin[1] < point[1] ≤ coordmax[1] &&
      coordmin[2] < point[2] ≤ coordmax[2] &&
      coordmin[3] < point[3] ≤ coordmax[3]
      return true
   else
      return false
   end
end

"""
    getcellinline(meta, point1, point2) -> cellids, distances, coords

Returns cell IDs, distances and coordinates for every cell in a line between two given
points `point1` and `point2`. May be improved later with preallocation!
"""
function getcellinline(meta::MetaVLSV, point1, point2)
   @unpack coordmin, coordmax, ncells = meta

   if !isInsideDomain(meta, point1)
      throw(DomainError(point1, "point location outside simulation domain!"))
   elseif !isInsideDomain(meta, point2)
      throw(DomainError(point2, "point location outside simulation domain!"))
   end

   cell_lengths = @inbounds [
      (coordmax[1]-coordmin[1])/ncells[1],
      (coordmax[2]-coordmin[2])/ncells[2],
      (coordmax[3]-coordmin[3])/ncells[3]]

   distances = [0.0]
   cellids = [getcell(meta, point1)]
   coords = point1
   ϵ = eps(Float32)
   unit_vector = @. (point2 - point1) / $hypot(point2 - point1 + ϵ...)
   p = point1

   while true
      cid = getcell(meta, p)
      amrlvl = getlevel(meta, cid)

      # Get the max and min cell boundaries
      min_bounds = getcellcoordinates(meta, cid) - 0.5*cell_lengths*0.5^amrlvl
      max_bounds = min_bounds + cell_lengths

      # Check which face we hit first
      coef_min = (min_bounds - p) ./ unit_vector
      coef_max = (max_bounds - p) ./ unit_vector

      # Negative coefficients indicates the opposite direction
      @inbounds for i = 1:3
         if unit_vector[i] == 0.0
            coef_min[i] = Inf
            coef_max[i] = Inf
         end
         if coef_min[i] ≤ 0  coef_min[i] = Inf end
         if coef_max[i] ≤ 0  coef_max[i] = Inf end
      end

      # Find the minimum distance from a boundary times a factor
      d = min(minimum(coef_min), minimum(coef_max)) * 1.00001

      coordnew = @. p + d*unit_vector

      dot(point2 - coordnew, unit_vector) ≥ 0 || break

      cellidnew = getcell(meta, coordnew)

      push!(cellids, cellidnew)
      coords = hcat(coords, coordnew)
      push!(distances, hypot(coordnew - point1...))

      p = coordnew
   end

   cellids, distances, coords
end

"""
    getslicecell(meta, sliceoffset;
       xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf) -> idlist, indexlist

Find the cell ids `idlist` which are needed to plot a 2d cut through of a 3d mesh, in a
direction given by non-Inf values for optional arguments at `sliceoffset`, which is always
within the range, and the `indexlist`, which is a mapping from original order to the cut
plane and can be used to select data onto the plane.
"""
function getslicecell(meta::MetaVLSV, sliceoffset;
   xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf)

   @unpack ncells, maxamr, cellid = meta

   if !isinf(xmin) && !isinf(xmax)
      minCoord = xmin; maxCoord = xmax; nsize = ncells[1]; idim = 1
   elseif !isinf(ymin) && !isinf(ymax)
      minCoord = ymin; maxCoord = ymax; nsize = ncells[2]; idim = 2
   elseif !isinf(zmin) && !isinf(zmax)
      minCoord = zmin; maxCoord = zmax; nsize = ncells[3]; idim = 3
   else
      throw(ArgumentError("Unspecified slice direction!"))
   end

   sliceratio = sliceoffset / (maxCoord - minCoord)
   @assert 0.0 ≤ sliceratio ≤ 1.0 "slice plane index out of bound!"

   # Find the ids
   nlen = 0
   ncell = prod(ncells)
   # number of cells up to each refinement level
   nStart = SVector{maxamr+2}(0, accumulate(+, (ncell*8^ilvl for ilvl = 0:maxamr))...)

   indexlist = Int[]
   idlist = Int[]

   @inbounds for ilvl = 0:maxamr
      nLow, nHigh = nStart[ilvl+1], nStart[ilvl+2]
      ids = cellid[nLow .< cellid .≤ nHigh]
      ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, ids)

      if idim == 1
         coords = ix
      elseif idim == 2
         coords = iy
      elseif idim == 3
         coords = iz
      end

      # Find the cut plane index for each refinement level (0-based)
      depth = floor(Int, sliceratio*nsize*2^ilvl)
      # Find the needed elements to create the cut and save the results
      elements = coords .== depth
      append!(indexlist, (nlen+1:nlen+length(ids))[elements])
      append!(idlist, ids[elements])

      nlen += length(ids)
   end

   idlist, indexlist
end

"""
    refineslice(meta, idlist, data, normal) -> Array

Generate scalar data on the finest refinement level given cellids `idlist` and variable
`data` on the slice perpendicular to `normal`.
"""
function refineslice(meta::MetaVLSV, idlist, data, normal)
   @unpack ncells, maxamr = meta

   if normal == :x
      dims = [ncells[2], ncells[3]] .* 2^maxamr
   elseif normal == :y
      dims = [ncells[1], ncells[3]] .* 2^maxamr
   elseif normal == :z
      dims = [ncells[1], ncells[2]] .* 2^maxamr
   end

   dpoints = zeros(eltype(data), dims...)

   # Create the plot grid
   ncell = prod(ncells)
   nHigh, nLow = ncell, 0

   @inbounds for i = 0:maxamr
      ids = idlist[nLow .< idlist .≤ nHigh]
      d = data[nLow .< idlist .≤ nHigh]

      ix, iy, iz = getindexes(i, ncells[1], ncells[2], nLow, ids)

      # Get the correct coordinate values and the widths for the plot
      if normal == :x
         a, b = iy, iz
      elseif normal == :y
         a, b = ix, iz
      elseif normal == :z
         a, b = ix, iy
      end

      # Insert the data values into dpoints
      refineRatio = 2^(maxamr - i)
      iRange = 0:refineRatio-1
      X, Y = ndgrid(iRange, iRange)

      coords = Array{Int64,3}(undef, 2, length(a), 2^(2*(maxamr-i)))
      @fastmath for ic in eachindex(a, b), ir = 1:2^(2*(maxamr-i))
         coords[1,ic,ir] = muladd(a[ic], refineRatio, 1+X[ir])
         coords[2,ic,ir] = muladd(b[ic], refineRatio, 1+Y[ir])
      end

      for ic in eachindex(d)
         dpoints[coords[1,ic,:],coords[2,ic,:]] .= d[ic]
      end

      nLow = nHigh
      nHigh += ncell*8^(i+1)
   end

   dpoints
end

"Compute every cell id's x, y and z indexes on the given refinement level (0-based)."
@inline function getindexes(ilevel, xcells, ycells, nCellUptoLowerLvl, ids)
   slicesize = xcells*ycells*4^ilevel

   iz = @. (ids - nCellUptoLowerLvl - 1) ÷ slicesize

   # number of ids up to the coordinate z in the refinement level ilevel
   idUpToZ = muladd.(iz, slicesize, nCellUptoLowerLvl)

   iy = @. (ids - idUpToZ - 1) ÷ (xcells*2^ilevel)
   ix = @. ids - idUpToZ - iy*xcells*2^ilevel - 1

   ix, iy, iz
end

@inline function getindexes(ilvl, xcells, ycells, nCellUptoLowerLvl, id::Int)
   slicesize = xcells*ycells*4^ilvl
   iz = (id - nCellUptoLowerLvl - 1) ÷ slicesize
   idUpToZ = muladd(iz, slicesize, nCellUptoLowerLvl)
   iy = (id - idUpToZ - 1) ÷ (xcells*2^ilvl)
   ix = id - idUpToZ - iy*xcells*2^ilvl - 1
   ix, iy, iz
end

"""
    getnearestcellwithvdf(meta, id) -> Int

Find the nearest spatial cell with VDF saved of a given cell `id` in the file `meta`.
"""
function getnearestcellwithvdf(meta::MetaVLSV, id)
   cells = getcellwithvdf(meta)
   isempty(cells) && throw(ArgumentError("No distribution saved in $(meta.name)"))
   coords = Matrix{Float32}(undef, 3, length(cells))
   @inbounds for i in eachindex(cells)
      coords[:,i] = getcellcoordinates(meta, cells[i])
   end
   coords_orig = getcellcoordinates(meta, id)
   d2 = sum((coords .- coords_orig).^2, dims=1)
   cells[argmin(d2)[2]]
end

"""
    getcellwithvdf(meta) -> cellids

Get all the cell IDs with VDF saved.
"""
getcellwithvdf(meta::MetaVLSV) =
   readmesh(meta.fid, meta.footer, "SpatialGrid", "CELLSWITHBLOCKS")


"Return the first cellid - 1 on `mylevel` given `ncells` on this level."
function get1stcell(mylevel, ncells)
   cid1st = 0
   for i = 0:mylevel-1
      cid1st += ncells*8^i
   end
   cid1st
end

fillmesh(meta::MetaVLSV, vars::AbstractString) = fillmesh(meta, [vars])

"""
    fillmesh(meta::MetaVLSV, vars; verbose=false) -> celldata, vtkGhostType

Fill the DCCRG mesh with quantity of `vars` on all refinement levels.
# Return arguments
- `celldata::Vector{Vector{Array}}`: data for each variable on each AMR level.
- `vtkGhostType::Array{UInt8}`: cell status (to be completed!).
"""
function fillmesh(meta::MetaVLSV, vars; verbose=false)
   @unpack cellid, maxamr, fid, footer, ncells = meta

   nvarvg = findall(!startswith("fg_"), vars)
   nv = length(vars)
   T = Vector{DataType}(undef, nv)
   offset = @MVector zeros(Int, nv)
   arraysize = @MVector zeros(Int, nv)
   dsize  = @MVector zeros(Int, nv)
   vsize  = @MVector zeros(Int, nv)
   @inbounds for i = 1:nv
      T[i], offset[i], arraysize[i], dsize[i], vsize[i] =
         getObjInfo(fid, footer, vars[i], "VARIABLE", "name")
   end

   Tout = copy(T)
   for i in eachindex(T)
      if T[i] == Float64 Tout[i] = Float32 end
   end

   @inbounds celldata =
      [[zeros(Tout[iv], vsize[iv], ncells[1]*2^i, ncells[2]*2^i, ncells[3]*2^i)
      for i = 0:maxamr] for iv in 1:nv]

   @inbounds vtkGhostType =
      [zeros(UInt8, ncells[1]*2^i, ncells[2]*2^i, ncells[3]*2^i) for i = 0:maxamr]

   if maxamr == 0
      @inbounds for iv = 1:nv
         celldata[iv][1][:] = readvariable(meta, vars[iv])
      end
      return celldata, vtkGhostType
   end

   # Find the ids
   ncell = prod(ncells)
   nLow, nHigh = 0, ncell
   cellidRaw = readvector(fid, footer, "CellID", "VARIABLE")

   @inbounds for ilvl = 0:maxamr
      verbose && @info "scanning AMR level $ilvl..."
      ids = cellid[nLow .< cellid .≤ nHigh]

      # indicate the condition of non-existing cells
      idrefined = setdiff(nLow+1:nHigh, ids)

      for id in idrefined
         ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, id) .+ 1
         vtkGhostType[ilvl+1][ix,iy,iz] = 8
      end

      rOffsetsRaw = indexin(ids, cellidRaw)

      if ilvl != maxamr
         for iv in nvarvg
            verbose && @info "reading variable $(vars[iv])..."
            a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv], offset[iv])
            dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
            data = @view dataRaw[:,rOffsetsRaw]

            fillcell!(iv, ilvl, ids, ncells, maxamr, nLow, celldata, data)
         end
      else # max amr level
         for (iv, var) = enumerate(vars)
            verbose && @info "reading variable $var..."
            if startswith(var, "fg_")
               celldata[iv][end][:] = readvariable(meta, var)
            else
               a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv],
                  offset[iv])
               dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
               data = @view dataRaw[:,rOffsetsRaw]

               fillcell!(iv, ids, ncells, maxamr, nLow, celldata, data)
            end
         end
      end
      nLow = nHigh
      nHigh += ncell*8^(ilvl+1)
   end

   celldata, vtkGhostType
end

function fillcell!(iv, ilvl, ids, ncells, maxamr, nLow, celldata, data)
   @inbounds for ilvlup = ilvl:maxamr
      r = 2^(ilvlup-ilvl) # ratio on refined level
      for c in eachindex(ids)
         ixr, iyr, izr = getindexes(ilvl, ncells[1], ncells[2], nLow, ids[c]) .* r
         for k = 1:r, j = 1:r, i = 1:r
            celldata[iv][ilvlup+1][:,ixr+i,iyr+j,izr+k] .= data[:,c]
         end
      end
   end
end

function fillcell!(iv, ids, ncells, maxamr, nLow, celldata, data)
   @inbounds for i in eachindex(ids)
      ix, iy, iz = getindexes(maxamr, ncells[1], ncells[2], nLow, ids[i]) .+ 1
      celldata[iv][end][:,ix,iy,iz] .= data[:,i]
   end
end


"""
    write_vtk(meta::MetaVLSV; kwargs...)
    write_vtk(filename; kwargs...)

Convert VLSV file to VTK format.
# Keyword arguments
- `vars=[""]`: select which variables to convert.
- `ascii=false`: output stored in ASCII or compressed binary format.
- `vti=false`: generate image files on the highest refinement level only.
- `verbose=false`: display logs during conversion.
"""
function write_vtk(meta::MetaVLSV; vars=[""], ascii=false, vti=false, verbose=false)
   @unpack ncells, maxamr, dcoord, coordmin = meta

   append = ascii ? false : true

   filedata = Vector{String}(undef, maxamr+1)
   @inbounds for i in 1:maxamr+1
      filedata[i] = meta.name[1:end-5]*"_$i.vti"
   end

   if isempty(vars[1])
      vars = meta.variable
      index_cellID = findfirst(==("CellID"), vars)
      if !isnothing(index_cellID) deleteat!(vars, index_cellID) end
   end

   data, vtkGhostType = fillmesh(meta, vars; verbose)

   if vti
      save_image(meta, meta.name[1:end-4]*"vti", vars, data, vtkGhostType[end], maxamr,
         append)
   else
      # Generate image file on each refinement level
      @inbounds for i in eachindex(vtkGhostType, filedata)
         fdata, ghost = filedata[i], vtkGhostType[i]
         save_image(meta, fdata, vars, data, ghost, i-1, append)
      end

      # Generate vthb file
      filemeta = meta.name[1:end-4]*"vthb"
      doc = XMLDocument()
      elm = ElementNode("VTKFile")
      setroot!(doc, elm)
      link!(elm, AttributeNode("type", "vtkOverlappingAMR"))
      link!(elm, AttributeNode("version", "1.1"))
      link!(elm, AttributeNode("byte_order", "LittleEndian")) # x86
      link!(elm, AttributeNode("header_type", "UInt64"))

      xamr = addelement!(elm, "vtkOverlappingAMR")

      origin = @sprintf "%f %f %f" coordmin...
      link!(xamr, AttributeNode("origin", origin))
      link!(xamr, AttributeNode("grid_description", "XYZ"))

      @inbounds for i = 0:maxamr
         xBlock = addelement!(xamr, "Block")
         link!(xBlock, AttributeNode("level", string(i)))
         spacing_str = @sprintf "%f %f %f" dcoord[1]/2^i dcoord[2]/2^i dcoord[3]/2^i
         link!(xBlock, AttributeNode("spacing", spacing_str))
         xDataSet = addelement!(xBlock, "DataSet")
         link!(xDataSet, AttributeNode("index", "0"))
         amr_box = [0, ncells[1]*2^i-1, 0, ncells[2]*2^i-1, 0, ncells[3]*2^i-1]
         box_str = @sprintf "%d %d %d %d %d %d" amr_box...
         link!(xDataSet, AttributeNode("amr_box", box_str))
         link!(xDataSet, AttributeNode("file", filedata[i+1]))
      end

      write(filemeta, doc)
   end

   return
end

write_vtk(filename; kwargs...) = write_vtk(load(filename); kwargs...)

"""
    save_image(meta::MetaVLSV, file, vars, data, vtkGhostType, level,
       ascii=false, append=true)

Save `data` of name `vars` at AMR `level` into VTK image file of name `file`.
# Arguments
- `file::String`: output file name.
- `vars::Vector{String}`: variable names to be saved.
- `data::Vector{Vector}`: data for all the variables on each refinement level.
- `vtkGhostType::Array{UInt8}`: array for visibility control.
- `level::Int`: refinement level (0-based).
- `ascii=false`: save output in ASCII or binary format.
- `append=true`: determines whether to append data at the end of file or do in-block writing.
"""
function save_image(meta::MetaVLSV, file, vars, data, vtkGhostType, level, ascii=false,
   append=true)
   @unpack coordmin, dcoord, ncells = meta
   origin = (coordmin[1], coordmin[2], coordmin[3])
   ratio = 2^level
   spacing = (dcoord[1] / ratio, dcoord[2] / ratio, dcoord[3] / ratio)

   vtk = vtk_grid(file, ncells[1]*ratio+1, ncells[2]*ratio+1, ncells[3]*ratio+1;
      origin, spacing, append, ascii)

   @inbounds for (iv, var) in enumerate(vars)
      vtk[var, VTKCellData()] = data[iv][level+1]
   end

   vtk["vtkGhostType", VTKCellData()] = vtkGhostType

   vtk_save(vtk)
end


"""
    issame(filename1, filename2, tol=1e-4) -> Bool

Check if two VLSV files are approximately identical.
"""
function issame(f1, f2, tol::AbstractFloat=1e-4; verbose=false)
   # 1st sanity check: minimal filesize difference
   if abs(filesize(f1) - filesize(f2)) / filesize(f2) > 1e-2
      verbose && println("The sizes of files are already quite different!")
      return false
   end

   meta1 = load(f1)
   meta2 = load(f2)
   varnames = meta1.variable
   strskip = r"CellID|rank|blocks"
   deleteat!(varnames, findall(endswith(strskip), varnames))

   isIdentical = true
   for vname in varnames
      v1 = readvariable(meta1, vname)
      v2 = readvariable(meta2, vname)

      s1, s2 = sum(v1), sum(v2)
      if abs(s1 - s2) > tol * abs(s1) && abs(s1 - s2) > tol * abs(s2)
         isIdentical = false
         verbose && println("$vname is quite different!")
         break
      end
   end
   verbose && isIdentical && println("$f1 and $f2 are identical under tolerance $tol.")
   close(meta1.fid)
   close(meta2.fid)
   return isIdentical
end