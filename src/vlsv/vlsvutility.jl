# Utility functions for processing VLSV data.

"""
    getcell(meta::MetaVLSV, location:::AbstractVector{<:AbstractFloat}) -> Int

Return cell ID containing the given spatial `location` in meter, excluding domain
boundaries. Only accept 3D location.
"""
function getcell(meta::MetaVLSV, loc::AbstractVector{<:AbstractFloat})
   (;coordmin, coordmax, dcoord, ncells, celldict, maxamr) = meta

   foreach( (i,comp) -> coordmin[i] < loc[i] < coordmax[i] ? nothing :
      error("$comp coordinate out of bound!"), 1:3, 'x':'z')

   indices = @inbounds ntuple(i -> round(Int, (loc[i] - coordmin[i]) ÷ dcoord[i]), Val(3))

   cid = @inbounds indices[1] + indices[2]*ncells[1] + indices[3]*ncells[1]*ncells[2] + 1

   ncells_lowerlevel = 0
   ncell = prod(ncells)

   @inbounds for ilevel = 0:maxamr
      haskey(celldict, cid) && break

      ncells_lowerlevel += (8^ilevel)*ncell
      ratio = 2^(ilevel+1)
      indices = ntuple(i -> floor(Int, (loc[i] - coordmin[i]) / dcoord[i] * ratio), Val(3))
      cid = ncells_lowerlevel + indices[1] +
         ratio*ncells[1]*indices[2] + ratio^2*ncells[1]*ncells[2]*indices[3] + 1
   end

   cid
end

"""
    getlevel(meta::MetaVLSV, cid::Int) -> Int

Return the AMR level of a given cell ID `cid`.
!!! warning
    This function does not check if the VLSV file of `meta` actually contains `cid`; it may
    be shadowed by refined children.
"""
function getlevel(meta::MetaVLSV, cid::Int)
   ncell = prod(meta.ncells)
   ilevel = 0
   c = Int(cid) - ncell
   while c > 0
      ilevel += 1
      c -= (8^ilevel)*ncell
   end

   ilevel
end

"""
    getparent(meta::MetaVLSV, cid::Int) -> Int

Return the parent cell ID of given child `cid`.
"""
function getparent(meta::MetaVLSV, cid::Int)
   @inbounds xcell, ycell = meta.ncells[1], meta.ncells[2]
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)
   parentlvl = mylvl - 1

   if parentlvl < 0
      throw(ArgumentError("Cell ID $cid has no parent cell!"))
   else
      # get the first cell ID on my level
      cid1st = get1stcell(mylvl, ncell)
      # get row and column sequence on my level (starting with 0)
      xcell <<= mylvl
      ycell <<= mylvl

      myseq = cid - cid1st
      ix = myseq % xcell
      iz = myseq ÷ (xcell*ycell)
      iy = (myseq - iz*xcell*ycell) ÷ xcell
      # indexes on the parent level
      ixparent = ix ÷ 2
      iyparent = iy ÷ 2
      izparent = iz ÷ 2

      # get the first cell ID on parent level
      cid1st -= ncell*8^parentlvl
      # get parent cell ID (may not exist!!!)
      parentid = cid1st + izparent*xcell*ycell÷4 + iyparent*xcell÷2 + ixparent
   end

   parentid
end

"""
    getchildren(meta::MetaVLSV, cid::Int) -> Vector{Int}

Return direct children of `cid`.
"""
function getchildren(meta::MetaVLSV, cid::Int)
   xcell, ycell, zcell = meta.ncells
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)

   # get the first cell ID on the my level
   cid1st = 1
   for i = 0:mylvl-1
      cid1st += ncell * 8^i
   end
   # get my row and column sequence on my level (starting with 0)
   xcell <<= mylvl
   ycell <<= mylvl

   myseq = cid - cid1st
   ix = myseq % xcell
   iz = myseq ÷ (xcell*ycell)
   iy = (myseq - iz*xcell*ycell) ÷ xcell

   # get the children sequences on the finer level
   ix *= 2
   iy *= 2
   iz *= 2

   nchildren = 2^ndims(meta)
   cid = Vector{Int}(undef, nchildren)
   # get the first cell ID on the finer level
   cid1st += ncell*8^mylvl
   ix_, iy_ = (ix, ix+1), (iy, iy+1)
   iz_ = zcell != 1 ? (iz, iz+1) : iz
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*xcell*ycell*4 + i[2]*xcell*2 + i[1]
   end

   (cid)
end

"""
    getsiblings(meta::MetaVLSV, cid::Int) -> Vector{Int}

Return sibling cells of a given `cid`, including itself.
"""
function getsiblings(meta::MetaVLSV, cid::Int)
   xcell, ycell, zcell = meta.ncells
   ncell = prod(meta.ncells)

   mylvl = getlevel(meta, cid)

   mylvl == 0 && throw(ArgumentError("CellID $cid is not a child cell!"))

   xcell = xcell << mylvl
   ycell = ycell << mylvl

   # 1st cell ID on my level
   cid1st = get1stcell(mylvl, ncell)

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
   cid = Vector{Int}(undef, nsiblings)
   ix_, iy_ = (ix, ix1), (iy, iy1)
   iz_ = zcell != 1 ? (iz, iz1) : iz
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*xcell*ycell + i[2]*xcell + i[1]
   end

   (cid)
end

"""
    isparent(meta::MetaVLSV, cid::Int) -> Bool

Check if `cid` is a parent cell.
"""
function isparent(meta::MetaVLSV, cid::Int)
   ncell_accum = get1stcell(meta.maxamr, prod(meta.ncells))

   !haskey(meta.celldict, cid) && 0 < cid < ncell_accum
end

"""
    getcellcoordinates(meta::MetaVLSV, cid::Int) -> SVector{3,Float64}

Return a given cell's spatial coordinates.
"""
function getcellcoordinates(meta::MetaVLSV, cid::Int)
   (;ncells, coordmin, coordmax) = meta
   cid -= 1 # for easy divisions

   ncells_refmax = collect(ncells)
   reflevel = 0
   subtraction = prod(ncells) * (2^reflevel)^3
   # sizes on the finest level
   while cid ≥ subtraction
      cid -= subtraction
      reflevel += 1
      subtraction *= 8
      ncells_refmax[1] *= 2
      ncells_refmax[2] *= 2
      ncells_refmax[3] *= 2
   end

   indices = @inbounds (
      cid % ncells_refmax[1],
      cid ÷ ncells_refmax[1] % ncells_refmax[2],
      cid ÷ (ncells_refmax[1] * ncells_refmax[2]) )

   coords = @inbounds @SVector [
      coordmin[i] + (indices[i] + 0.5) * (coordmax[i] - coordmin[i]) / ncells_refmax[i]
      for i = 1:3]

   coords
end

"""
    getvcellcoordinates(meta::MetaVLSV, vcellids::Vector; species="proton")

Return velocity cells' coordinates of `species` and `vcellids`.
"""
function getvcellcoordinates(meta::MetaVLSV, vcellids::Vector{Int32};
   species::String="proton")
   (;vblocks, vblock_size, dv, vmin) = meta.meshes[species]

   bsize = prod(vblock_size)
   blockid = @. vcellids ÷ bsize
   # Get block coordinates
   blockInd = [(
      bid % vblocks[1],
      bid ÷ vblocks[1] % vblocks[2],
      bid ÷ (vblocks[1] * vblocks[2]) )
      for bid in blockid]
   blockCoord = [(
      bInd[1] * dv[1] * vblock_size[1] + vmin[1],
      bInd[2] * dv[2] * vblock_size[2] + vmin[2],
      bInd[3] * dv[3] * vblock_size[3] + vmin[3] )
      for bInd in blockInd]

   # Get cell indices
   vcellblockids = @. vcellids % bsize
   cellidxyz = [(
      cid % vblock_size[1],
      cid ÷ vblock_size[1] % vblock_size[2],
      cid ÷ (vblock_size[1] * vblock_size[2]) )
      for cid in vcellblockids]

   # Get cell coordinates
   cellCoords = [SVector(0.0f0, 0.0f0, 0.0f0) for _ in vcellblockids]
   @inbounds @simd for i in eachindex(vcellblockids)
      cellCoords[i] = [blockCoord[i][j] + (cellidxyz[i][j] + 0.5) * dv[j] for j in 1:3]
   end

   cellCoords
end

"""
    getdensity(meta, VDF; species="proton")
    getdensity(meta, vcellf; species="proton")
    getdensity(vmesh::VMeshInfo, vcellf)

Get density from `VDF` of `species` associated with `meta`, n = ∫ f(r,v) dV. Alternatively,
one can directly pass `vcellids` as original indices of nonzero VDFs and `vcellf` as their
corresponding values.
"""
function getdensity(meta::MetaVLSV, VDF::Array{T};
   species::String="proton") where T <: AbstractFloat

   (;dv) = meta.meshes[species]

   n = sum(VDF) * convert(T, prod(dv))
end

function getdensity(vmesh::VMeshInfo, vcellf::Vector{T}) where T <: AbstractFloat
   n = sum(vcellf) * convert(T, prod(vmesh.dv))
end

getdensity(meta::MetaVLSV, vcellf::Vector{<:AbstractFloat}; species::String="proton") =
   getdensity(meta.meshes[species], vcellf)

"""
    getvelocity(meta, VDF; species="proton")
    getvelocity(meta, vcellids, vcellf; species="proton")
    getvelocity(vmesh::VMeshInfo, vcellids, vcellf)

Get bulk velocity from `VDF` of `species`, u = ∫ v * f(r,v) dV / n. Alternatively, one can
directly pass `vcellids`, `vcellf`, as in [`getdensity`](@ref).
"""
function getvelocity(meta::MetaVLSV, VDF::Array{T};
   species::String="proton") where T <: AbstractFloat

   (;dv, vmin) = meta.meshes[species]
   u = @SVector zeros(T, 3)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]
      u += VDF[i,j,k] .* SVector(vx, vy, vz)
   end

   n = sum(VDF)

   u /= n

   u
end

function getvelocity(vmesh::VMeshInfo, vcellids::Vector{Int32}, vcellf::Vector{T}) where
   T <: AbstractFloat

   (;vblock_size, vblocks, dv, vmin) = vmesh
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   slicez = vsize[1]*vsize[2]
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   sliceCz = vblock_size[1]*vblock_size[2]

   u = @SVector zeros(T, 3)

   @inbounds @simd for ic in eachindex(vcellids)
      id = findindex(vcellids[ic], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
      i = id % vsize[1]
      j = id % slicez ÷ vsize[1]
      k = id ÷ slicez

      vx = vmin[1] + (i + 0.5f0)*dv[1]
      vy = vmin[2] + (j + 0.5f0)*dv[2]
      vz = vmin[3] + (k + 0.5f0)*dv[3]
      u += vcellf[ic] .* SVector(vx, vy, vz)
   end

   n = sum(vcellf)

   u /= n

   u
end

getvelocity(meta::MetaVLSV, vcellids::Vector{Int32}, vcellf::Vector{T};
   species::String="proton") where T <: AbstractFloat =
   getvelocity(meta.meshes[species], vcellids, vcellf)

"""
    getpressure(meta, VDF; species="proton")
    getpressure(meta, vcellids, vcellf; species="proton")
    getpressure(vmesh::VMeshInfo, vcellids, vcellf)

Get pressure tensor (6 components: Pxx, Pyy, Pzz, Pyz, Pzx, Pxy) of `species` from `VDF`
associated with `meta`, pᵢⱼ = m/3 * ∫ (v - u)ᵢ(v - u)ⱼ * f(r,v) dV.
Alternatively, one can directly pass `vcellids`, `vcellf`, as in [`getdensity`](@ref).
"""
function getpressure(meta::MetaVLSV, VDF::Array{T};
   species::String="proton") where T <: AbstractFloat

   (;dv, vmin) = meta.meshes[species]
   p = @SVector zeros(T, 6)

   u = getvelocity(meta, VDF; species)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]

      p += VDF[i,j,k] .* SVector(
         (vx - u[1])*(vx - u[1]),
         (vy - u[2])*(vy - u[2]),
         (vz - u[3])*(vz - u[3]),
         (vy - u[2])*(vz - u[3]),
         (vx - u[1])*(vz - u[3]),
         (vx - u[1])*(vy - u[2]))
   end

   factor = mᵢ * convert(T, prod(dv))
   p *= factor

   p
end

function getpressure(vmesh::VMeshInfo, vcellids::Vector{Int32}, vcellf::Vector{T}) where
   T <: AbstractFloat

   (;vblock_size, vblocks, dv, vmin) = vmesh
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   slicez = vsize[1]*vsize[2]
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   sliceCz = vblock_size[1]*vblock_size[2]

   u = getvelocity(vmesh, vcellids, vcellf)

   p = @SVector zeros(T, 6)

   @inbounds @simd for ic in eachindex(vcellids)
      id = findindex(vcellids[ic], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
      i = id % vsize[1]
      j = id % slicez ÷ vsize[1]
      k = id ÷ slicez

      vx = vmin[1] + (i + 0.5f0)*dv[1]
      vy = vmin[2] + (j + 0.5f0)*dv[2]
      vz = vmin[3] + (k + 0.5f0)*dv[3]
      p += vcellf[ic] .* SVector(
         (vx - u[1])*(vx - u[1]),
         (vy - u[2])*(vy - u[2]),
         (vz - u[3])*(vz - u[3]),
         (vy - u[2])*(vz - u[3]),
         (vx - u[1])*(vz - u[3]),
         (vx - u[1])*(vy - u[2]))
   end

   factor = mᵢ * convert(T, prod(dv))
   p *= factor

   p
end

getpressure(meta::MetaVLSV, vcellids::Vector{Int32}, vcellf::Vector{<:AbstractFloat};
   species::String="proton") =
   getpressure(meta.meshes[species], vcellids, vcellf)

"""
   getheatfluxvector(meta, VDF; species="proton")
   getheatfluxvector(meta, vcellids, vcellf; species="proton")
   getheatfluxvector(vmesh::VMeshInfo, vcellids, vcellf)

Get heat flux vector (3 components) of `species` from `VDF` associated with `meta`,
qᵢ = m/2 * ∫ (v - u)²(v - u)ᵢ * f(r,v) dV. Alternatively, one can directly pass `vcellids`,
`vcellf`, as in [`getdensity`](@ref).
"""
function getheatfluxvector(meta::MetaVLSV, VDF::Array{T}; species::String="proton") where
   T <: AbstractFloat

   (;dv, vmin) = meta.meshes[species]
   q = @SVector zeros(T, 3)

   u = getvelocity(meta, VDF; species)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]

      q += VDF[i,j,k] .* SVector(
         (vx - u[1])*(vx - u[1])*(vx - u[1]),
         (vy - u[2])*(vy - u[2])*(vy - u[2]),
         (vz - u[3])*(vz - u[3])*(vz - u[3]))
   end

   factor = mᵢ * convert(T, prod(dv))
   q *= factor

   q
end

function getheatfluxvector(vmesh::VMeshInfo, vcellids::Vector{Int32}, vcellf::Vector{T}
   ) where T <: AbstractFloat

   (;vblock_size, vblocks, dv, vmin) = vmesh
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   slicez = vsize[1]*vsize[2]
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   sliceCz = vblock_size[1]*vblock_size[2]

   u = getvelocity(vmesh, vcellids, vcellf)

   q = @SVector zeros(T, 3)

   @inbounds @simd for ic in eachindex(vcellids)
      id = findindex(vcellids[ic], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
      i = id % vsize[1]
      j = id % slicez ÷ vsize[1]
      k = id ÷ slicez

      vx = vmin[1] + (i + 0.5f0)*dv[1]
      vy = vmin[2] + (j + 0.5f0)*dv[2]
      vz = vmin[3] + (k + 0.5f0)*dv[3]
      q += vcellf[ic] .* SVector(
         (vx - u[1])*(vx - u[1])*(vx - u[1]),
         (vy - u[2])*(vy - u[2])*(vy - u[2]),
         (vz - u[3])*(vz - u[3])*(vz - u[3]))
  end

   factor = mᵢ * convert(T, prod(dv))
   q *= factor

   q
end

getheatfluxvector(meta::MetaVLSV, vcellids::Vector{Int32}, vcellf::Vector{<:AbstractFloat};
   species::String="proton") =
   getheatfluxvector(meta.meshes[species], vcellids, vcellf)

"Get the original vcell index without blocks from raw vcell index `i` (0-based)."
@inline function findindex(i::Int32, vblocks::NTuple{3, Int}, vblock_size::NTuple{3,Int},
   blocksize::Int, vsize::NTuple{3, Int}, sliceBz::Int, sliceCz::Int)
   iB = i ÷ blocksize
   iBx = iB % vblocks[1]
   iBy = iB % sliceBz ÷ vblocks[1]
   iBz = iB ÷ sliceBz
   iCellInBlock = i % blocksize
   iCx = iCellInBlock % vblock_size[1]
   iCy = iCellInBlock % sliceCz ÷ vblock_size[1]
   iCz = iCellInBlock ÷ sliceCz
   iBCx = iBx*vblock_size[1] + iCx
   iBCy = iBy*vblock_size[2] + iCy
   iBCz = iBz*vblock_size[3] + iCz

   iOrigin = iBCz*vsize[1]*vsize[2] + iBCy*vsize[1] + iBCx
end

"""
    reorder(vmesh::VMeshInfo, vcellids::Vector{Int32}) -> vcellids_origin

Reorder vblock-organized VDF indexes into x-->y-->z indexes. `vcellids` are raw indices
of nonzero VDFs ordered by blocks.
"""
function reorder(vmesh::VMeshInfo, vcellids::Vector{Int32})
   (;vblock_size, vblocks) = vmesh
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   sliceCz = vblock_size[1]*vblock_size[2]

   vcellids_origin = similar(vcellids)
   # IDs are 0-based
   @inbounds @simd for i in eachindex(vcellids)
      vcellids_origin[i] = 1 +
         findindex(vcellids[i], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
   end

   vcellids_origin
end

"""
    reconstruct(vmesh::VMeshInfo, vcellids::Vector{Int32}, vcellf::Vector{<:AbstractFloat})

Reconstruct the full VDFs in 3D. `vcellids` are raw indices of nonzero VDFs ordered by
blocks, and `vcellf` are the corresponding values in each cell.
"""
function reconstruct(vmesh::VMeshInfo, vcellids::Vector{Int32},
   vcellf::Vector{<:AbstractFloat})
   (;vblock_size, vblocks) = vmesh
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   sliceCz = vblock_size[1]*vblock_size[2]
   # Reconstruct the full velocity space
   VDF = zeros(Float32, vsize)
   # Raw IDs are 0-based
   @inbounds @simd for i in eachindex(vcellids)
      j = 1 +
         findindex(vcellids[i], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
      VDF[j] = vcellf[i]
   end

   VDF
end

"""
    getmaxwellianity(meta, VDF; species="proton")
    getmaxwellianity(meta, vcellids, vcellf; species="proton")

Obtain the Maxwellian similarity factor -log(1/(2n) * ∫ |f - g| dv), where `f` is the VDF
from Vlasiator and `g` is the analytical Maxwellian distribution that generates the same
density as `f`. The value ranges from [0, +∞], with 0 meaning not Maxwellian-distributed at
all, and +∞ a perfect Maxwellian distribution.
Alternatively, one can pass original `vcellids` and `vcellf` directly.
"""
function getmaxwellianity(meta::MetaVLSV, VDF::Array{<:AbstractFloat};
   species::String="proton")
   (;dv, vmin) = meta.meshes[species]

   n = getdensity(meta, VDF)
   u = getvelocity(meta, VDF)
   P = getpressure(meta, VDF)
   p = (P[1] + P[2] + P[3]) / 3
   T = p / (n *kB) # temperature from scalar pressure

   ϵₘ = zero(eltype(VDF))
   vth2⁻¹ = mᵢ / (2kB*T)

   @inbounds for k in axes(VDF,3), j in axes(VDF,2), i in axes(VDF,1)
      vx = vmin[1] + (i - 0.5f0)*dv[1]
      vy = vmin[2] + (j - 0.5f0)*dv[2]
      vz = vmin[3] + (k - 0.5f0)*dv[3]
      dv2 = (vx - u[1])^2 + (vy - u[2])^2 + (vz - u[3])^2

      g = n * √(vth2⁻¹/π) * (vth2⁻¹/π) * ℯ^(-vth2⁻¹*dv2)
      ϵₘ += abs(VDF[i,j,k] - g)
   end

   ϵₘ = -log(0.5 / n * convert(eltype(VDF), prod(dv)) * ϵₘ)
end

function getmaxwellianity(meta::MetaVLSV, vcellids::Vector{Int32},
   vcellf::Vector{<:AbstractFloat}; species::String="proton")

   (;vblock_size, vblocks, dv, vmin) = meta.meshes[species]
   vsize = @inbounds ntuple(i -> vblock_size[i] * vblocks[i], Val(3))
   slicez = vsize[1]*vsize[2]
   blocksize = prod(vblock_size)
   sliceBz = vblocks[1]*vblocks[2]
   sliceCz = vblock_size[1]*vblock_size[2]

   n = getdensity(meta, vcellf)
   u = getvelocity(meta, vcellids, vcellf)
   P = getpressure(meta, vcellids, vcellf)
   p = (P[1] + P[2] + P[3]) / 3
   T = p / (n *kB) # temperature from scalar pressure

   ϵₘ = zero(eltype(vcellf))
   vth2⁻¹ = mᵢ / (2kB*T)

   @inbounds @simd for ic in eachindex(vcellids)
      id = findindex(vcellids[ic], vblocks, vblock_size, blocksize, vsize, sliceBz, sliceCz)
      i = id % vsize[1]
      j = id % slicez ÷ vsize[1]
      k = id ÷ slicez

      vx = vmin[1] + (i + 0.5f0)*dv[1]
      vy = vmin[2] + (j + 0.5f0)*dv[2]
      vz = vmin[3] + (k + 0.5f0)*dv[3]
      dv2 = (vx - u[1])^2 + (vy - u[2])^2 + (vz - u[3])^2

      g = n * √(vth2⁻¹/π) * (vth2⁻¹/π) * ℯ^(-vth2⁻¹*dv2)
      ϵₘ += abs(vcellf[ic] - g)
   end

   ϵₘ = -log(0.5 / n * convert(eltype(vcellf), prod(dv)) * ϵₘ)
end

function isInsideDomain(meta::MetaVLSV, point::Vector{<:Real})
   (;coordmin, coordmax) = meta

   if coordmin[1] < point[1] ≤ coordmax[1] &&
      coordmin[2] < point[2] ≤ coordmax[2] &&
      coordmin[3] < point[3] ≤ coordmax[3]
      return true
   else
      return false
   end
end

"""
    getcellinline(meta, point1::Vector, point2::Vector) -> cellids, distances, coords

Returns cell IDs, distances and coordinates for every cell in a line between two given
points `point1` and `point2`. TODO: preallocation?
"""
function getcellinline(meta::MetaVLSV, point1::Vector{T}, point2::Vector{T}) where T
   (;coordmin, coordmax, ncells) = meta

   if !isInsideDomain(meta, point1)
      throw(DomainError(point1, "point location outside simulation domain!"))
   elseif !isInsideDomain(meta, point2)
      throw(DomainError(point2, "point location outside simulation domain!"))
   end

   dcell = @inbounds ntuple(i -> (coordmax[i] - coordmin[i]) / ncells[i], Val(3))

   distances = [zero(T)]
   cellids = [getcell(meta, point1)]
   coords = [SVector{3}(point1)]
   ϵ = eps(T)
   unit_vector = @. (point2 - point1) / $norm(point2 - point1 + ϵ)
   p = coords[1]
   coef_min = Vector{T}(undef, 3)
   coef_max = Vector{T}(undef, 3)

   @inbounds while true
      cid = getcell(meta, p)
      amrlvl = getlevel(meta, cid)

      # Get the max and min cell boundaries
      Δ = dcell.*0.5.^amrlvl
      min_bounds = getcellcoordinates(meta, cid) .- 0.5.*Δ
      max_bounds = min_bounds .+ Δ

      # Check which face we hit first
      @. coef_min = (min_bounds - p) / unit_vector
      @. coef_max = (max_bounds - p) / unit_vector

      # Negative coefficients indicates the opposite direction
      for i = 1:3
         if unit_vector[i] == 0.0
            coef_min[i] = Inf
            coef_max[i] = Inf
         end
         if coef_min[i] ≤ 0  coef_min[i] = Inf end
         if coef_max[i] ≤ 0  coef_max[i] = Inf end
      end

      # Find the minimum distance from a boundary times a factor
      d = min(minimum(coef_min), minimum(coef_max)) * 1.00001

      coordnew = SVector(
         p[1] + d*unit_vector[1],
         p[2] + d*unit_vector[2],
         p[3] + d*unit_vector[3])

      dot(point2 .- coordnew, unit_vector) ≥ 0 || break

      cellidnew = getcell(meta, coordnew)

      push!(cellids, cellidnew)
      push!(coords, coordnew)
      push!(distances, norm(coordnew .- point1))

      p = coordnew
   end

   cellids, distances, coords
end

"""
    getslicecell(meta, sliceoffset, dir, minCoord, maxCoord) -> idlist, indexlist

Find the cell IDs `idlist` which are needed to plot a 2d cut through of a 3d mesh, in a
direction `dir` at `sliceoffset`, and the `indexlist`, which is a mapping from original
order to the cut plane and can be used to select data onto the plane.
"""
function getslicecell(meta::MetaVLSV, sliceoffset::Float64, dir::Int,
   minCoord::Float64, maxCoord::Float64)
   dir ∉ (1,2,3) && @error "Unknown slice direction $dir"
   (;ncells, maxamr, celldict) = meta

   nsize = ncells[dir]
   sliceratio = sliceoffset / (maxCoord - minCoord)
   0.0 ≤ sliceratio ≤ 1.0 || error("slice plane index out of bound!")

   # Find the ids
   nlen = 0
   ncell = prod(ncells)
   # number of cells up to each refinement level
   nStart = (vcat(0, accumulate(+, (ncell*8^ilvl for ilvl = 0:maxamr))))

   indexlist = Int[]
   idlist = Int[]

   cellidsorted = sort(collect(keys(celldict)))

   @inbounds for ilvl = 0:maxamr
      nLow, nHigh = nStart[ilvl+1], nStart[ilvl+2]
      idfirst_ = searchsortedfirst(cellidsorted, nLow+1)
      idlast_  = searchsortedlast(cellidsorted, nHigh)

      ids = @view cellidsorted[idfirst_:idlast_]

      ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, ids)

      coords =
         if dir == 1
            ix
         elseif dir == 2
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
    refineslice(meta, idlist::Vector{Int}, data::AbstractVector, normal::Symbol) -> Vector
    refineslice(meta, idlist::Vector{Int}, data::AbstractMatrix, normal::Symbol) -> Matrix
    refineslice(meta, idlist::Vector{Int}, data::AbstractArray, dir::Int)

Generate data on the finest refinement level given cellids `idlist` and variable `data` on
the slice perpendicular to `normal`. If `data` is 2D, then the first dimension is treated as
vector components.
"""
function refineslice(meta::MetaVLSV, idlist::Vector{Int}, data::AbstractVector,
   normal::Symbol)
   (;ncells, maxamr) = meta

   dims = _getdim2d(ncells, maxamr, normal)

   dpoints = zeros(eltype(data), dims...)

   # Create the plot grid
   ncell = prod(ncells)
   nHigh, nLow = ncell, 0

   @inbounds for i = 0:maxamr
      idfirst_ = searchsortedfirst(idlist, nLow+1)
      idlast_  = searchsortedlast(idlist, nHigh)

      ids = @view idlist[idfirst_:idlast_]
      d = @view data[idfirst_:idlast_]

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

      coords = [(0, 0) for _ in a, _ in 1:2^(2*(maxamr-i))]

      @inbounds for ir = 1:2^(2*(maxamr-i)), ic in eachindex(a, b)
         @fastmath coords[ic,ir] = (muladd(a[ic], refineRatio, 1+X[ir]),
                                    muladd(b[ic], refineRatio, 1+Y[ir]) )
      end

      for ir = 1:2^(2*(maxamr-i)), ic in eachindex(d)
         dpoints[ coords[ic,ir]... ] = d[ic]
      end

      nLow = nHigh
      nHigh += ncell*8^(i+1)
   end

   dpoints
end

function refineslice(meta::MetaVLSV, idlist::Vector{Int}, data::AbstractArray, dir::Int)
   if dir == 1
      refineslice(meta, idlist, data, :x)
   elseif dir == 2
      refineslice(meta, idlist, data, :y)
   elseif dir == 3
      refineslice(meta, idlist, data, :z)
   else
      error("Normal direction index $normal out of range!")
   end
end

function refineslice(meta::MetaVLSV, idlist::Vector{Int}, data::AbstractMatrix,
   normal::Symbol)
   dims = _getdim2d(meta.ncells, meta.maxamr, normal)

   dout = zeros(eltype(data), size(data,1), dims...)

   for i in axes(data,1)
      dout[i,:,:] = @views refineslice(meta, idlist, data[i,:], normal)
   end

   dout
end

@inline function _getdim2d(ncells::NTuple{3, Int}, maxamr::Int, normal::Symbol)
   ratio = 2^maxamr
   if normal == :x
      i1, i2 = 2, 3
   elseif normal == :y
      i1, i2 = 1, 3
   elseif normal == :z
      i1, i2 = 1, 2
   end
   dims = (ncells[i1]*ratio, ncells[i2]*ratio)
end

"Compute x, y and z indexes of all cell `ids` on the given refinement level (0-based)."
@inline function getindexes(ilevel::Int, xcells::Int, ycells::Int, nCellUptoLowerLvl::Int,
   ids::AbstractVector{Int})

   ratio = 2^ilevel
   slicesize = xcells*ycells*ratio^2

   iz = @. (ids - nCellUptoLowerLvl - 1) ÷ slicesize
   iy = similar(iz)
   ix = similar(iz)
   @inbounds for i in eachindex(ids, iz)
      # number of ids up to the coordinate z in the refinement level ilevel
      idUpToZ = muladd(iz[i], slicesize, nCellUptoLowerLvl)
      iy[i] = (ids[i] - idUpToZ - 1) ÷ (xcells*ratio)
      ix[i] = ids[i] - idUpToZ - iy[i]*xcells*ratio - 1
   end

   ix, iy, iz
end

"Compute x, y and z index of cell `id` on a given refinement level `ilevel`(0-based)."
@inline function getindexes(ilevel::Int, xcells::Int, ycells::Int, nCellUptoLowerLvl::Int,
   id::Int)

   ratio = 2^ilevel
   slicesize = xcells*ycells*ratio^2
   iz = (id - nCellUptoLowerLvl - 1) ÷ slicesize
   idUpToZ = muladd(iz, slicesize, nCellUptoLowerLvl)
   iy = (id - idUpToZ - 1) ÷ (xcells*ratio)
   ix = id - idUpToZ - iy*xcells*ratio - 1

   ix, iy, iz
end

"""
    getnearestcellwithvdf(meta, id::Int) -> Int

Find the nearest spatial cell with VDF saved of a given cell `id` associated with `meta`.
"""
function getnearestcellwithvdf(meta::MetaVLSV, id::Int)
   cells = getcellwithvdf(meta)
   isempty(cells) && throw(ArgumentError("No distribution saved in $(meta.name)"))
   coords = [SVector(0.0f0, 0.0f0, 0.0f0) for _ in cells]
   @inbounds for i in eachindex(cells)
      coords[i] = getcellcoordinates(meta, cells[i])
   end
   coords_orig = getcellcoordinates(meta, id)
   d2 = [sum((c .- coords_orig).^2) for c in coords]

   cells[argmin(d2)]
end

"""
    getcellwithvdf(meta, species::String="proton") -> cellids

Get all the cell IDs with VDF saved associated with `meta`.
"""
function getcellwithvdf(meta::MetaVLSV, species::String="proton")
   (; fid, nodeVLSV) = meta

   local cellsWithVDF, nblock_C

   for node in nodeVLSV.cellwithVDF
      if node["name"] == species
         asize = Parsers.parse(Int, node["arraysize"])
         cellsWithVDF = Vector{Int}(undef, asize)
         offset = Parsers.parse(Int, nodecontent(node))
         seek(fid, offset)
         read!(fid, cellsWithVDF)
         break
      end
   end

   for node in nodeVLSV.cellblocks
      if node["name"] == species
         asize = Parsers.parse(Int, node["arraysize"])
         dsize = Parsers.parse(Int, node["datasize"])
         offset = Parsers.parse(Int, nodecontent(node))
         nblock_C = dsize == 4 ?
            Vector{Int32}(undef, asize) : Vector{Int}(undef, asize)
         seek(fid, offset)
         read!(fid, nblock_C)
         break
      end
   end

   innerBCCells = findall(==(0), nblock_C)

   deleteat!(cellsWithVDF, innerBCCells)

   cellsWithVDF
end

"Return the first cell ID on `mylevel` given `ncells` on this level."
get1stcell(mylevel::Int, ncells::Int) = ncells * (8^mylevel - 1) ÷ 7 + 1

## Downsampling

"""
    downsample_fg(meta::MetaVLSV, v_fg::Array)
    downsample_fg(meta::MetaVLSV, var::String)

Downsample a field solver array `v_fg` to the spatial grid associated with `meta`.
"""
function downsample_fg(meta::MetaVLSV, v_fg::Array)
   v_vg = zeros(eltype(v_fg), (3, length(meta.cellindex)))

   for (cid, i) in meta.celldict
      v_vg[:,i] = downsample_fg_cell(meta, v_fg, cid)
   end

   v_vg[:,meta.cellindex]
end

downsample_fg(meta::MetaVLSV, var::String) = downsample_fg(meta, meta[var])

"Return a field solver grid subarray contained inside spatial cell `cid`."
function downsample_fg_cell(meta::MetaVLSV, v_fg::Array, cid::Int)
   v_sub_fg = get_fg_array_cell(meta, v_fg, cid)

   dropdims(mean(v_sub_fg, dims=(2,3,4)), dims=(2,3,4))
end

"Return the field solver grid cell indexes containing `coords` (low-inclusive)."
function get_fg_indices(meta::MetaVLSV, coords::SVector{3, Float64})
   dx = meta.dcoord ./ 2^meta.maxamr
   ri = @. Int(((coords - meta.coordmin) ÷ dx) + 1)
   sz = meta.ncells .* 2 .^meta.maxamr
   if any(i -> i < 1, ri) || any(i -> i[1] > i[2], zip(ri, sz))
      error("fsgrid index out of bounds!")
   end

   ri
end

"""
    get_fg_array_cell(meta::MetaVLSV, v_fg::Array, cid::Int)

Return a subarray of the field solver grid array, corresponding to the fsgrid covered by
the spatial cell ID `cid`.
"""
function get_fg_array_cell(meta::MetaVLSV, v_fg::Array, cid::Int)
   il, ih = get_fg_indices_cell(meta, cid)
   v_fg[:,il[1]:ih[1],il[2]:ih[2],il[3]:ih[3]]
end

"Returns a slice tuple of fsgrid indices that are contained in the spatial cell `cid`."
function get_fg_indices_cell(meta::MetaVLSV, cid::Int)
   dx = meta.dcoord ./ 2^(getlevel(meta, cid) + 1)
   mid = getcellcoordinates(meta, cid)

   il, ih = get_fg_indices_subvolume(meta, mid .- dx, mid .+ dx)
end

"""
    get_fg_indices_subvolume(meta::MetaVLSV, lower, upper, tol::Float64=1e-3)

Get indices for subarray of fsgrid variables, in a cuboid defined by `lower` and `upper`
vertices. This is used for mapping a set of fsgrid cells to a given DCCRG cell.
Shift the corners (`lower`, `upper`) inward by a distance controlled by `tol`. If direct
low-inclusive behaviour is required, `tol` shall be set to 0.
"""
function get_fg_indices_subvolume(meta::MetaVLSV, lower::SVector{3, Float64},
   upper::SVector{3, Float64}, tol::Float64=1e-3)
   ϵ = @. meta.dcoord / 2^meta.maxamr * tol
   il = get_fg_indices(meta, lower .+ ϵ)
   iu = get_fg_indices(meta, upper .- ϵ)

   il, iu
end

## Upsampling

"""
    read_variable_as_fg(meta::MetaVLSV, var::String)

Interpolate DCCRG variable `var` to field solver grid size.
This is an alternative method to [`fillmesh`](@ref), but not optimized for performance.
"""
function read_variable_as_fg(meta::MetaVLSV, var::String)
   sz = meta.ncells .* 2 .^meta.maxamr
   data = readvariable(meta, var, false)
   if eltype(data) == Float64; data = Float32.(data); end
   cellid = Vlasiator.getcellid(meta.fid, meta.nodeVLSV.var)
   if ndims(data) == 2
      v_fg = zeros(eltype(data), (size(data,1), sz[1], sz[2], sz[3]))
      for (i, cid) in enumerate(cellid)
         upsample_fsgrid_subarray!(meta, data[:,i], cid, v_fg)
      end
   else
      v_fg = zeros(eltype(data), (sz[1], sz[2], sz[3]))
      for (i, cid) in enumerate(cellid)
         upsample_fsgrid_subarray!(meta, data[i], cid, v_fg)
      end
   end

   v_fg
end

"Set the elements of the fsgrid array to the value of corresponding cell ID `cid`."
function upsample_fsgrid_subarray!(meta::MetaVLSV, data, cid::Int, v_fg::Array{T, 4}
   ) where T 
   il, ih = get_fg_indices_cell(meta, cid)
   v_fg[:,il[1]:ih[1],il[2]:ih[2],il[3]:ih[3]] .= data
end

function upsample_fsgrid_subarray!(meta::MetaVLSV, data, cid::Int, v_fg::Array{T, 3}
   ) where T 
   il, ih = get_fg_indices_cell(meta, cid)
   v_fg[il[1]:ih[1],il[2]:ih[2],il[3]:ih[3]] .= data
end

fillmesh(meta::MetaVLSV, vars::String) = fillmesh(meta, [vars])

"""
    fillmesh(meta::MetaVLSV, vars::Vector{String};
       skipghosttype=true, maxamronly=false, verbose=false) -> celldata, vtkGhostType

Fill the DCCRG mesh with quantity of `vars` on all refinement levels.
# Return arguments
- `celldata::Vector{Vector{Array}}`: data for each variable on each AMR level.
- `vtkGhostType::Array{UInt8}`: cell status (to be completed!).
"""
function fillmesh(meta::MetaVLSV, vars::Vector{String};
   skipghosttype::Bool=true, maxamronly::Bool=false, verbose::Bool=false)

   (;maxamr, fid, nodeVLSV, ncells, celldict) = meta

   nvarvg = findall(!startswith("fg_"), vars)
   nv = length(vars)
   T = Vector{DataType}(undef, nv)
   offset = Vector{Int}(undef, nv)
   arraysize = similar(offset)
   vsize = similar(offset)
   @inbounds for i = 1:nv
      T[i], offset[i], arraysize[i], _, vsize[i] = getvarinfo(nodeVLSV.var, vars[i])
   end

   Tout = copy(T)
   for i in eachindex(T)
      if T[i] == Float64 Tout[i] = Float32 end
   end

   @inbounds celldata = if !maxamronly
      [[zeros(Tout[iv], vsize[iv], ncells[1] << i, ncells[2] << i, ncells[3] << i)
      for i = 0:maxamr] for iv in 1:nv]
   else
      [[zeros(Tout[iv], vsize[iv],
      ncells[1] << maxamr, ncells[2] << maxamr, ncells[3] << maxamr)] for iv in 1:nv]
   end

   @inbounds vtkGhostType = if !skipghosttype
      [zeros(UInt8, ncells[1] << i, ncells[2] << i, ncells[3] << i) for i = 0:maxamr]
   else
      [UInt8[]]
   end

   if maxamr == 0
      @inbounds for iv = 1:nv
         celldata[iv][1][:] = readvariable(meta, vars[iv])
      end
      return celldata, vtkGhostType
   end

   # Find the ids
   ncell = prod(ncells)
   nLow, nHigh = 0, ncell
   cellidsorted = sort(collect(keys(celldict)))

   @inbounds for ilvl = 0:maxamr
      verbose && @info "scanning AMR level $ilvl..."

      idfirst_ = searchsortedfirst(cellidsorted, nLow+1)
      idlast_  = searchsortedlast(cellidsorted, nHigh)

      ids = @view cellidsorted[idfirst_:idlast_]

      if !skipghosttype
         # Mark non-existing cells due to refinement
         @simd for id in nLow+1:nHigh
            if isempty(searchsorted(ids, id))
               ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, id) .+ 1
               vtkGhostType[ilvl+1][ix,iy,iz] = 8
            end
         end
      end

      rOffsetsRaw = [celldict[id] for id in ids]

      ix, iy, iz = getindexes(ilvl, ncells[1], ncells[2], nLow, ids)

      if ilvl != maxamr
         for iv in nvarvg
            verbose && @info "reading variable $(vars[iv])..."
            a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv], offset[iv])
            dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
            data = @view dataRaw[:,rOffsetsRaw]
            if !maxamronly
               _fillcell!(ids, ix, iy, iz, celldata[iv], data, ilvl, maxamr)
            else
               _fillcell!(ids, ix, iy, iz, celldata[iv][end], data, ilvl, maxamr)
            end
         end
      else # max refinement level
         for (iv, var) in enumerate(vars)
            verbose && @info "reading variable $var..."
            if startswith(var, "fg_")
               celldata[iv][end][:] = readvariable(meta, var)
            else
               a = mmap(fid, Vector{UInt8}, sizeof(T[iv])*vsize[iv]*arraysize[iv],
                  offset[iv])
               dataRaw = reshape(reinterpret(T[iv], a), vsize[iv], arraysize[iv])
               data = @view dataRaw[:,rOffsetsRaw]

               _fillcell!(ids, ix, iy, iz, celldata[iv][end], data)
            end
         end
      end
      nLow = nHigh
      nHigh += ncell*8^(ilvl+1)
   end

   celldata, vtkGhostType
end

function _fillcell!(ids::AbstractVector{Int}, ix::Vector{Int}, iy::Vector{Int},
   iz::Vector{Int}, dataout::Vector, datain::AbstractArray, ilvl::Int, maxamr::Int)
   @inbounds for ilvlup = ilvl:maxamr
      r = 2^(ilvlup-ilvl) # ratio on refined level
      for c in eachindex(ids)
         for k = 1:r, j = 1:r, i = 1:r
            _fillcelldata!(dataout[ilvlup+1], datain, ix[c]*r+i, iy[c]*r+j, iz[c]*r+k, c)
         end
      end
   end
end

function _fillcell!(ids::AbstractVector{Int}, ix::Vector{Int}, iy::Vector{Int},
   iz::Vector{Int}, dataout::Array, datain::AbstractArray, ilvl::Int, maxamr::Int)
   r = 2^(maxamr-ilvl) # ratio on refined level
   @inbounds for c in eachindex(ids)
      for k = 1:r, j = 1:r, i = 1:r
         _fillcelldata!(dataout, datain, ix[c]*r+i, iy[c]*r+j, iz[c]*r+k, c)
      end
   end
end

function _fillcell!(ids::AbstractVector{Int}, ix::Vector{Int}, iy::Vector{Int},
   iz::Vector{Int}, dataout::Array, datain::AbstractArray)
   @inbounds for c in eachindex(ids)
      _fillcelldata!(dataout, datain, ix[c]+1, iy[c]+1, iz[c]+1, c)
   end
end

@inline function _fillcelldata!(dataout::Array, datain::AbstractArray,
   i::Int, j::Int, k::Int, index::Int)
   @inbounds @simd for icomp in axes(datain, 1)
      dataout[icomp,i,j,k] = datain[icomp,index]
   end
end

"""
    write_vtk(meta::MetaVLSV; kwargs...)
    write_vtk(file::AbstractString; kwargs...)

Convert VLSV file to VTK format.
# Keywords
- `vars::Vector{String}=[""]`: select which variables to convert.
- `ascii::Bool=false`: output stored in ASCII or compressed binary format.
- `maxamronly::Bool=false`: generate image files on the highest refinement level only.
- `box::Vector`: selected box range in 3D.
- `outdir::String=""`: output directory.
- `verbose::Bool=false`: display logs during conversion.
"""
function write_vtk(meta::MetaVLSV; vars::Vector{String}=[""], ascii::Bool=false,
   maxamronly::Bool=false, skipghosttype::Bool=false, verbose::Bool=false,
   box::Vector{Float64}=[-Inf, Inf, -Inf, Inf, -Inf, Inf], outdir::AbstractString="")

   (;ncells, maxamr, dcoord, coordmin) = meta

   append = ascii ? false : true

   filedata = Vector{String}(undef, maxamr+1)
   @inbounds for i in 1:maxamr+1
      filedata[i] = outdir*meta.name[1:end-5]*"_$i.vti"
   end

   if isempty(vars[1])
      vars = meta.variable
      cellid_ = findfirst(==("CellID"), vars)
      if !isnothing(cellid_) deleteat!(vars, cellid_) end
   end

   data, vtkGhostType = fillmesh(meta, vars; skipghosttype, maxamronly, verbose)

   if maxamronly
      save_image(meta, outdir*meta.name[1:end-4]*"vti", vars, data, vtkGhostType[end],
         maxamr, append; box)
   else
      # Generate image file on each refinement level
      @inbounds for i in eachindex(vtkGhostType, filedata)
         fdata, ghost = filedata[i], vtkGhostType[i]
         save_image(meta, fdata, vars, data, ghost, i-1, append)
      end

      # Generate vthb file
      filemeta = outdir*meta.name[1:end-4]*"vthb"
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
         amr_box = (0, ncells[1]*2^i-1, 0, ncells[2]*2^i-1, 0, ncells[3]*2^i-1)
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
       ascii=false, append=true, box=[-Inf, Inf, -Inf, Inf, -Inf, Inf])

Save `data` of name `vars` at AMR `level` into VTK image file of name `file`.
# Arguments
- `file::String`: output file name.
- `vars::Vector{String}`: variable names to be saved.
- `data::Vector{Vector{Array}}`: data for all the variables on each refinement level.
- `vtkGhostType::Array{UInt8}`: array for visibility control.
- `level::Int`: refinement level (0-based).
- `ascii::Bool=false`: save output in ASCII or binary format.
- `append::Bool=true`: determines whether to append data at the end of file or do in-block
writing.
- `box::Vector`: selected box range in 3D.
"""
function save_image(meta::MetaVLSV, file::String, vars::Vector{String},
   data, vtkGhostType::Array{UInt8}, level::Int, ascii::Bool=false, append::Bool=true;
   box::Vector{<:AbstractFloat}=[-Inf, Inf, -Inf, Inf, -Inf, Inf])

   (;coordmin, coordmax, dcoord, ncells) = meta
   ratio = 2^level
   spacing = (dcoord[1] / ratio, dcoord[2] / ratio, dcoord[3] / ratio)
   # Only max amr level is stored
   if length(data[1]) == 1
      level = 0
   end

   if all(isinf.(box)) # full domain
      origin = (coordmin[1], coordmin[2], coordmin[3])

      vtk = vtk_grid(file, ncells[1]*ratio+1, ncells[2]*ratio+1, ncells[3]*ratio+1;
         origin, spacing, append, ascii)

      @inbounds for (iv, var) in enumerate(vars)
         vtk[var, VTKCellData()] = data[iv][level+1]
      end

      vtk["vtkGhostType", VTKCellData()] = vtkGhostType
   else # selected box region
      xrange = box[1]:spacing[1]:box[2]
      yrange = box[3]:spacing[2]:box[4]
      zrange = box[5]:spacing[3]:box[6]
      vtk = vtk_grid(file, xrange, yrange, zrange; append, ascii)

      if box[1] < coordmin[1] || box[3] < coordmin[2] || box[5] < coordmin[3] ||
         box[2] > coordmax[1] || box[4] > coordmax[2] || box[6] > coordmax[3]
         throw(ArgumentError("Selected box range $box out of bound!"))
      end

      imin = Int((box[1] - coordmin[1]) ÷ spacing[1]) + 1
      imax = imin + length(xrange) - 2
      jmin = Int((box[3] - coordmin[2]) ÷ spacing[2]) + 1
      jmax = jmin + length(yrange) - 2
      kmin = Int((box[5] - coordmin[3]) ÷ spacing[3]) + 1
      kmax = kmin + length(zrange) - 2

      @inbounds for (iv, var) in enumerate(vars)
         vtk[var, VTKCellData()] = data[iv][level+1][:,imin:imax,jmin:jmax,kmin:kmax]
      end
      vtk["vtkGhostType", VTKCellData()] = vtkGhostType[imin:imax,jmin:jmax,kmin:kmax]
   end

   vtk_save(vtk)
end

"""
    write_vlsv(filein, fileout, newvars::Vector{Tuple{Vector, String, VarInfo}};
       force=false)

Generate a new VLSV `fileout` based on `filein`, with `newvars` added.
`force=true` overwrites the existing `fileout`.
"""
function write_vlsv(filein::AbstractString, fileout::AbstractString,
   newvars::Vector{Tuple{VecOrMat, String, VarInfo}}; force::Bool=false)

   if isfile(fileout) && !force
      error("Output target $fileout exists!")
   end

   fid = open(filein)
   endian_offset = 8 # First 8 bytes indicate big-endian or else
   seek(fid, endian_offset)
   # Obtain the offset of the XML footer
   offset = read(fid, Int)
   # Store all non-footer part as raw data
   raw_data = Vector{UInt8}(undef, offset)

   seekstart(fid)
   readbytes!(fid, raw_data, offset)
   # Read input VLSV file footer
   doc = read(fid, String) |> parsexml
   footer = doc |> root
   close(fid)
   # Get new variables' offsets
   offsets = accumulate(+,
      [offset, [sizeof(newvars[i][1]) for i in eachindex(newvars)[1:end-1]]...])
   # Create new children for footer
   for i in eachindex(newvars, offsets)
      elm = addelement!(footer, "VARIABLE", string(offsets[i]))

      a1 = AttributeNode("arraysize", string(length(newvars[i][1])))
      a2 = AttributeNode("datasize", string(sizeof(eltype(newvars[i][1]))))
      a3 =
         if eltype(newvars[i][1]) <: Signed
            AttributeNode("datatype", "int")
         elseif eltype(newvars[i][1]) <: AbstractFloat
            AttributeNode("datatype", "float")
         elseif eltype(newvars[i][1]) <: Unsigned
            AttributeNode("datatype", "uint")
         end
      a4 = AttributeNode("mesh", "SpatialGrid")
      a5 = AttributeNode("name", newvars[i][2])
      a6 = AttributeNode("unit", newvars[i][3].unit)
      a7 = AttributeNode("unitConversion", newvars[i][3].unitConversion)
      a8 = AttributeNode("unitLaTeX", newvars[i][3].unitLaTeX)
      a9 = AttributeNode("variableLaTeX", newvars[i][3].variableLaTeX)
      a10 =
         if ndims(newvars[i][1]) == 1
            AttributeNode("vectorsize", "1")
         else
            AttributeNode("vectorsize", string(size(newvars[i][1], 1)))
         end

      for attributenode in (a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)
         link!(elm, attributenode)
      end
   end
   # Write to fileout
   open(fileout, "w") do io
      write(io, @view raw_data[1:8]) # endianness
      # Compute footer offset
      totalnewsize = 0
      for var in newvars
         totalnewsize += sizeof(var[1])
      end
      write(io, offset+totalnewsize) # record new footer offset
      write(io, @view raw_data[17:end]) # copy original data
      for var in newvars
         write(io, var[1])
      end
      write(io, string(footer), '\n')
   end

   return
end

"""
    issame(file1, file2, tol=1e-4; strict=true, verbose=false) -> Bool

Check if two VLSV files `file1` and `file2` are approximately identical, under relative
tolerance `tol`. If `strict=true`, the file size difference should be less than 1%.
"""
function issame(f1::AbstractString, f2::AbstractString, tol::AbstractFloat=1e-4;
   strict::Bool=true, verbose::Bool=false)
   # 1st sanity check: minimal filesize difference
   if strict && abs(filesize(f1) - filesize(f2)) / filesize(f2) > 1e-2
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
      end
   end
   verbose && isIdentical && println("$f1 and $f2 are identical under tolerance $tol.")
   close(meta1.fid)
   close(meta2.fid)

   return isIdentical
end