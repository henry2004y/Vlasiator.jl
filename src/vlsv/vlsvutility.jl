# Utility functions for processing VLSV data.

"""
    getcell(meta, location) -> UInt

Return cell ID containing the given spatial `location`, excluding domain boundaries.
Only accept 3D location.
"""
function getcell(meta::MetaVLSV, loc)
   @unpack coordmin, coordmax, dcoord, ncells, cellid, maxamr = meta

   coordmin[1] < loc[1] < coordmax[1] || error("x coordinate out of bound!")
   coordmin[2] < loc[2] < coordmax[2] || error("y coordinate out of bound!")
   coordmin[3] < loc[3] < coordmax[3] || error("z coordinate out of bound!")

   dx, dy, dz = dcoord

   # Get cell indices
   indices = @inbounds SVector(
      round(UInt, (loc[1] - coordmin[1]) ÷ dx),
      round(UInt, (loc[2] - coordmin[2]) ÷ dy),
      round(UInt, (loc[3] - coordmin[3]) ÷ dz) )
   # Get cell id
   cid = @inbounds indices[1] + indices[2]*ncells[1] + indices[3]*ncells[1]*ncells[2] + 1

   ncells_lowerlevel = UInt(0)
   ncell = prod(ncells)

   @inbounds for ilevel = 0:maxamr
      cid in cellid && break

      ncells_lowerlevel += 2^(3*ilevel)*ncell

      ratio = 2^(ilevel+1)

      indices = SVector(
         floor(UInt, (loc[1] - coordmin[1]) / dx * ratio),
         floor(UInt, (loc[2] - coordmin[2]) / dy * ratio),
         floor(UInt, (loc[3] - coordmin[3]) / dz * ratio) )

      cid = ncells_lowerlevel + indices[1] +
         ratio*ncells[1]*indices[2] + ratio^2*ncells[1]*ncells[2]*indices[3] + 1
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
   c = Int(cid) - ncell
   while c > 0
      ilevel += 1
      c -= 2^(3*ilevel)*ncell
   end
   ilevel
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
   ix_, iy_ = (ix, ix+1), (iy, iy+1)
   iz_ = zcell != 1 ? (iz, iz+1) : iz
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
   ix_, iy_ = (ix, ix1), (iy, iy1)
   iz_ = zcell != 1 ? (iz, iz1) : iz
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
    getcellcoordinates(meta, cid) -> SVector{Float64}

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

   indices = SVector(
      cid % xcell,
      cid ÷ xcell % ycell,
      cid ÷ (xcell*ycell) )

   coords = @inbounds SVector(
      coordmin[1] + (indices[1] + 0.5) * (coordmax[1] - coordmin[1]) / xcell,
      coordmin[2] + (indices[2] + 0.5) * (coordmax[2] - coordmin[2]) / ycell,
      coordmin[3] + (indices[3] + 0.5) * (coordmax[3] - coordmin[3]) / zcell )

   coords
end

"""
    getvcellcoordinates(meta, vcellids, species="proton")

Return velocity cells' coordinates of `species` and `vcellids`.
"""
function getvcellcoordinates(meta::MetaVLSV, vcellids; species="proton")
   @unpack vblocks, vblock_size, dv, vmin = meta.meshes[species]

   bsize = prod(vblock_size)
   blockid = @. vcellids ÷ bsize
   # Get block coordinates
   blockInd = [SVector(
      bid % vblocks[1],
      bid ÷ vblocks[1] % vblocks[2],
      bid ÷ (vblocks[1] * vblocks[2]) )
      for bid in blockid]
   blockCoord = [SVector(
      bInd[1] * dv[1] * vblock_size[1] + vmin[1],
      bInd[2] * dv[2] * vblock_size[2] + vmin[2],
      bInd[3] * dv[3] * vblock_size[3] + vmin[3] )
      for bInd in blockInd]

   # Get cell indices
   vcellblockids = @. vcellids % bsize
   cellidxyz = [SVector(
      cid % vblock_size[1],
      cid ÷ vblock_size[1] % vblock_size[2],
      cid ÷ (vblock_size[1] * vblock_size[2]) )
      for cid in vcellblockids]

   # Get cell coordinates
   cellCoords = [zeros(SVector{3, Float32}) for _ in vcellblockids]
   @inbounds for i in eachindex(vcellblockids)
      cellCoords[i] = SVector(blockCoord[i][1] + (cellidxyz[i][1] + 0.5) * dv[1],
                              blockCoord[i][2] + (cellidxyz[i][2] + 0.5) * dv[2],
                              blockCoord[i][3] + (cellidxyz[i][3] + 0.5) * dv[3] )
   end
   cellCoords
end

"""
    getdensity(meta, VDF; species="proton")
    getdensity(meta, vcellids, vcellf; species="proton")

Get density from VDF, n = ∫ f(r,v) dV.
"""
function getdensity(meta::MetaVLSV, VDF; species="proton")
   @unpack dv = meta.meshes[species]
   n = zero(eltype(VDF))

   @inbounds @simd for f in VDF
      n += f
   end
   n * convert(eltype(VDF), prod(dv))
end

function getdensity(meta::MetaVLSV, vcellids, vcellf; species="proton")
   @unpack dv = meta.meshes[species]

   n = zero(eltype(vcellf))

   @inbounds @simd for f in vcellf
      n += f
   end
   n * convert(eltype(vcellf), prod(dv))
end

"""
    getvelocity(meta, VDF; species="proton")
    getvelocity(meta, vcellids, vcellf; species="proton")

Get bulk velocity from VDF, u = ∫ v * f(r,v) dV / n.
"""
function getvelocity(meta::MetaVLSV, VDF; species="proton")
   @unpack dv, vmin = meta.meshes[species]
   u = zeros(eltype(VDF), 3)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]
      u[1] += vx*VDF[i,j,k]
      u[2] += vy*VDF[i,j,k]
      u[3] += vz*VDF[i,j,k]
   end

   n = zero(eltype(VDF))
   @inbounds @simd for f in VDF
      n += f
   end

   u ./= n
end

function getvelocity(meta::MetaVLSV, vcellids, vcellf; species="proton")
   @unpack dv, vmin = meta.meshes[species]

   VDF = flatten(meta.meshes[species], vcellids, vcellf)

   u = getvelocity(meta, VDF; species)
end

"""
    getpressure(VDF)

Get pressure tensor from VDF, pᵢⱼ = m/3 * ∫ (v - u)ᵢ(v - u)ⱼ * f(r,v) dV.
"""
function getpressure(meta::MetaVLSV, VDF; species="proton")
   @unpack dv, vmin = meta.meshes[species]
   p = zeros(eltype(VDF), 6)

   u = getvelocity(meta, VDF; species)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]

      p[1] += (vx - u[1])*(vx - u[1])*VDF[i,j,k]
      p[2] += (vy - u[2])*(vy - u[2])*VDF[i,j,k]
      p[3] += (vz - u[3])*(vz - u[3])*VDF[i,j,k]
      p[4] += (vy - u[2])*(vz - u[3])*VDF[i,j,k]
      p[5] += (vx - u[1])*(vz - u[3])*VDF[i,j,k]
      p[6] += (vx - u[1])*(vy - u[2])*VDF[i,j,k]
   end

   p .*= mᵢ / 3 * convert(eltype(VDF), prod(dv))
end

"""
    flatten(vmesh::VMeshInfo, vcellids, vcellf)

Flatten vblock-organized VDFs into x-->y-->z ordered 3D VDFs. `vcellids` are local indices
of nonzero VDFs and `vcellf` are their corresponding values.
"""
function flatten(vmesh::VMeshInfo, vcellids, vcellf)
   vblock_size, vblocks = vmesh.vblock_size, vmesh.vblocks
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   vsizex, vsizey = vblock_size[1] * vblocks[1], vblock_size[2] * vblocks[2]
   sliceCz = vblock_size[1]*vblock_size[2]
   # Reconstruct the full velocity space
   VDFraw = zeros(Float32,
      vmesh.vblock_size[1]*vmesh.vblocks[1],
      vmesh.vblock_size[2]*vmesh.vblocks[2],
      vmesh.vblock_size[3]*vmesh.vblocks[3])
   # Fill nonzero values
   @inbounds for i in eachindex(vcellids)
      VDFraw[vcellids[i]+1] = vcellf[i]
   end

   VDF = zeros(Float32,
      vmesh.vblock_size[1]*vmesh.vblocks[1],
      vmesh.vblock_size[2]*vmesh.vblocks[2],
      vmesh.vblock_size[3]*vmesh.vblocks[3])

   @inbounds for i in eachindex(VDF)
      iB = (i - 1) ÷ blocksize
      iBx = iB % vblocks[1]
      iBy = iB % sliceBz ÷ vblocks[1]
      iBz = iB ÷ sliceBz
      iCellInBlock = (i - 1) % blocksize
      iCx = iCellInBlock % vblock_size[1]
      iCy = iCellInBlock % sliceCz ÷ vblock_size[1]
      iCz = iCellInBlock ÷ sliceCz
      iBCx = iBx*vblock_size[1] + iCx
      iBCy = iBy*vblock_size[2] + iCy
      iBCz = iBz*vblock_size[3] + iCz
      iF = iBCz*vsizex*vsizey + iBCy*vsizex + iBCx + 1
      VDF[iF] = VDFraw[i]
   end
   VDF
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

   cell_lengths = @inbounds SVector(
      (coordmax[1] - coordmin[1]) / ncells[1],
      (coordmax[2] - coordmin[2]) / ncells[2],
      (coordmax[3] - coordmin[3]) / ncells[3] )

   distances = [0.0]
   cellids = [getcell(meta, point1)]
   coords = point1
   ϵ = eps(Float32)
   unit_vector = @. (point2 - point1) / $norm(point2 - point1 + ϵ)
   p = point1
   coef_min = MVector(0.0, 0.0, 0.0)
   coef_max = MVector(0.0, 0.0, 0.0)

   while true
      cid = getcell(meta, p)
      amrlvl = getlevel(meta, cid)

      # Get the max and min cell boundaries
      min_bounds = getcellcoordinates(meta, cid) - 0.5*cell_lengths*0.5^amrlvl
      max_bounds = min_bounds + cell_lengths

      # Check which face we hit first
      @. coef_min = (min_bounds - p) / unit_vector
      @. coef_max = (max_bounds - p) / unit_vector

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

      coordnew = @inbounds SVector(p[1] + d*unit_vector[1],
         p[2] + d*unit_vector[2], p[3] + d*unit_vector[3] )

      dot(point2 - coordnew, unit_vector) ≥ 0 || break

      cellidnew = getcell(meta, coordnew)

      push!(cellids, cellidnew)
      coords = hcat(coords, coordnew)
      push!(distances, norm(coordnew .- point1))

      p = coordnew
   end

   cellids, distances, coords
end

"""
    getslicecell(meta, sliceoffset, idim, minCoord, maxCoord) -> idlist, indexlist

Find the cell ids `idlist` which are needed to plot a 2d cut through of a 3d mesh, in a
direction `idim` at `sliceoffset`, and the `indexlist`, which is a mapping from original
order to the cut plane and can be used to select data onto the plane.
"""
function getslicecell(meta::MetaVLSV, sliceoffset, idim, minCoord, maxCoord)
   idim ∉ (1,2,3) && @error "Unknown slice direction $idim"
   @unpack ncells, maxamr, cellid, cellindex = meta

   nsize = ncells[idim]
   sliceratio = sliceoffset / (maxCoord - minCoord)
   0.0 ≤ sliceratio ≤ 1.0 || error("slice plane index out of bound!")

   # Find the ids
   nlen = 0
   ncell = prod(ncells)
   # number of cells up to each refinement level
   nStart = SVector{maxamr+2}(vcat(0, accumulate(+, (ncell*8^ilvl for ilvl = 0:maxamr))))

   indexlist = Int[]
   idlist = UInt[]

   cellidsorted = cellid[cellindex]

   @inbounds for ilvl = 0:maxamr
      nLow, nHigh = nStart[ilvl+1], nStart[ilvl+2]
      ids = cellidsorted[nLow .< cellidsorted .≤ nHigh]
      ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, ids)

      coords =
         if idim == 1
            ix
         elseif idim == 2
            iy
         else # 3
            iz
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

   dims = let ratio = 2^maxamr
      if normal == :x
         i1, i2 = 2, 3
      elseif normal == :y
         i1, i2 = 1, 3
      elseif normal == :z
         i1, i2 = 1, 2
      end
      SVector(ncells[i1]*ratio, ncells[i2]*ratio)
   end

   dpoints = zeros(eltype(data), dims[1], dims[2])

   # Create the plot grid
   ncell = prod(ncells)
   nHigh, nLow = ncell, 0

   @inbounds for i = 0:maxamr
      ids = idlist[nLow .< idlist .≤ nHigh]
      d = data[nLow .< idlist .≤ nHigh]

      ix, iy, iz = getindexes(i, ncells[1], ncells[2], nLow, ids)

      # Get the correct coordinate values and the widths for the plot
      a, b =
         if normal == :x
            iy, iz
         elseif normal == :y
            ix, iz
         elseif normal == :z
            ix, iy
         end

      # Insert the data values into dpoints
      refineRatio = 2^(maxamr - i)
      iRange = 0:refineRatio-1
      X, Y = ndgrid(iRange, iRange)

      coords = [SVector(0, 0) for _ in a, _ in 1:2^(2*(maxamr-i))]
      @fastmath for ir = 1:2^(2*(maxamr-i)), ic in eachindex(a, b)
         coords[ic,ir] = SVector(muladd(a[ic], refineRatio, 1+X[ir]),
                                 muladd(b[ic], refineRatio, 1+Y[ir]) )
      end

      for ic in eachindex(d), ir = 1:2^(2*(maxamr-i))
         dpoints[ coords[ic,ir][1], coords[ic,ir][2] ] = d[ic]
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
   iy = similar(iz)
   ix = similar(iz)
   @inbounds for i in eachindex(ids, iz)
      # number of ids up to the coordinate z in the refinement level ilevel
      idUpToZ = muladd(iz[i], slicesize, nCellUptoLowerLvl)
      iy[i] = (ids[i] - idUpToZ - 1) ÷ (xcells*2^ilevel)
      ix[i] = ids[i] - idUpToZ - iy[i]*xcells*2^ilevel - 1
   end
   ix, iy, iz
end

@inline function getindexes(ilvl, xcells, ycells, nCellUptoLowerLvl, id::Integer)
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
   coords = [zeros(SVector{3, Float32}) for _ in cells]
   @inbounds for i in eachindex(cells)
      coords[i] = getcellcoordinates(meta, cells[i])
   end
   coords_orig = getcellcoordinates(meta, id)
   d2 = [sum((c - coords_orig).^2) for c in coords]
   cells[argmin(d2)]
end

"""
    getcellwithvdf(meta) -> cellids

Get all the cell IDs with VDF saved.
"""
function getcellwithvdf(meta::MetaVLSV)
   fid, footer = meta.fid, meta.footer
   cellsWithVDF = readmesh(fid, footer, "SpatialGrid", "CELLSWITHBLOCKS")::Vector{UInt}
   nblock_C = readmesh(fid, footer, "SpatialGrid", "BLOCKSPERCELL")::Vector{UInt32}

   innerBCCells = findall(==(0), nblock_C)

   deleteat!(cellsWithVDF, innerBCCells)
   cellsWithVDF
end


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
   @unpack maxamr, fid, footer, ncells, cellid, cellindex = meta

   nvarvg = findall(!startswith("fg_"), vars)
   nv = length(vars)
   T = Vector{DataType}(undef, nv)
   offset = @MVector zeros(Int, nv)
   arraysize = @MVector zeros(Int, nv)
   vsize  = @MVector zeros(Int, nv)
   @inbounds for i = 1:nv
      T[i], offset[i], arraysize[i], _, vsize[i] =
         getObjInfo(footer, vars[i], "VARIABLE", "name")
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
   cellidsorted = cellid[cellindex]

   @inbounds for ilvl = 0:maxamr
      verbose && @info "scanning AMR level $ilvl..."
      ids = cellidsorted[nLow .< cellidsorted .≤ nHigh]

      # indicate the condition of non-existing cells
      idrefined = setdiff(nLow+1:nHigh, ids)

      for id in idrefined
         ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, id) .+ 1
         vtkGhostType[ilvl+1][ix,iy,iz] = 8
      end

      rOffsetsRaw = indexin(ids, cellid)

      if ilvl != maxamr
         for iv in nvarvg
            verbose && @info "reading variable $(vars[iv])..."
            a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv], offset[iv])
            dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
            data = @view dataRaw[:,rOffsetsRaw]

            fillcell!(ilvl, ids, ncells, maxamr, nLow, celldata[iv], data)
         end
      else # max refinement level
         for (iv, var) = enumerate(vars)
            verbose && @info "reading variable $var..."
            if startswith(var, "fg_")
               celldata[iv][end][:] = readvariable(meta, var)
            else
               a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv],
                  offset[iv])
               dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
               data = @view dataRaw[:,rOffsetsRaw]

               fillcell!(ids, ncells, maxamr, nLow, celldata[iv][end], data)
            end
         end
      end
      nLow = nHigh
      nHigh += ncell*8^(ilvl+1)
   end

   celldata, vtkGhostType
end

function fillcell!(ilvl, ids, ncells, maxamr, nLow, dataout, datain)
   @inbounds for ilvlup = ilvl:maxamr
      r = 2^(ilvlup-ilvl) # ratio on refined level
      for c in eachindex(ids)
         ixr, iyr, izr = getindexes(ilvl, ncells[1], ncells[2], nLow, ids[c]) .* r
         for k = 1:r, j = 1:r, i = 1:r
            _fillcelldata!(dataout[ilvlup+1], datain, ixr+i, iyr+j, izr+k, c)
         end
      end
   end
end

function fillcell!(ids, ncells, maxamr, nLow, dataout, datain)
   @inbounds for i in eachindex(ids)
      ix, iy, iz = getindexes(maxamr, ncells[1], ncells[2], nLow, ids[i]) .+ 1
      _fillcelldata!(dataout, datain, ix, iy, iz, i)
   end
end

@inline function _fillcelldata!(dataout, datain, i, j, k, index)
   @inbounds for icomp in axes(datain,1)
      dataout[icomp,i,j,k] = datain[icomp,index]
   end
end

"""
    write_vtk(meta::MetaVLSV; kwargs...)
    write_vtk(file; kwargs...)

Convert VLSV file to VTK format.
# Keyword arguments
- `vars=[""]`: select which variables to convert.
- `ascii=false`: output stored in ASCII or compressed binary format.
- `maxamronly=false`: generate image files on the highest refinement level only.
- `verbose=false`: display logs during conversion.
"""
function write_vtk(meta::MetaVLSV; vars=[""], ascii=false, maxamronly=false, verbose=false)
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

   if maxamronly
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

      origin = @sprintf "%f %f %f" coordmin[1] coordmin[2] coordmin[3]
      link!(xamr, AttributeNode("origin", origin))
      link!(xamr, AttributeNode("grid_description", "XYZ"))

      @inbounds for i = 0:maxamr
         xBlock = addelement!(xamr, "Block")
         link!(xBlock, AttributeNode("level", string(i)))
         spacing_str = @sprintf "%f %f %f" dcoord[1]/2^i dcoord[2]/2^i dcoord[3]/2^i
         link!(xBlock, AttributeNode("spacing", spacing_str))
         xDataSet = addelement!(xBlock, "DataSet")
         link!(xDataSet, AttributeNode("index", "0"))
         amr_box = SA[0, ncells[1]*2^i-1, 0, ncells[2]*2^i-1, 0, ncells[3]*2^i-1]
         box_str = @sprintf("%d %d %d %d %d %d", amr_box[1], amr_box[2], amr_box[3],
            amr_box[4], amr_box[5], amr_box[6])
         link!(xDataSet, AttributeNode("amr_box", box_str))
         link!(xDataSet, AttributeNode("file", filedata[i+1]))
      end

      write(filemeta, doc)
   end

   return
end

write_vtk(file; kwargs...) = write_vtk(load(file); kwargs...)

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
    issame(file1, file2, tol=1e-4; verbose=false) -> Bool

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