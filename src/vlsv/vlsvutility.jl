# Utility functions for processing VLSV data.

using LinearAlgebra: dot
using WriteVTK, Printf

const qₑ = -1.60217662e-19  # electron charge, [C]
const mₑ = 9.10938356e-31   # electron mass, [kg]
const qᵢ = 1.60217662e-19   # proton mass, [C]
const mᵢ = 1.673557546e-27  # proton mass, [kg]
const c  = 3e8              # speed of light, [m/s]
const μ₀ = 4π*1e-7          # Vacuum permeability, [H/m]
const kB = 1.38064852e-23   # Boltzmann constant, [m²kg/(s²K)] 
const Re = 6.371e6          # Earth radius, [m]

export getcell, getslicecell, getlevel, getmaxamr, refineslice, getcellcoordinates,
   getchildren, getparent, isparent, getsiblings,
   getcellinline, getnearestcellwithvdf, write_vtk, compare

"""
    getcell(meta, location) -> Int

Return cell ID containing the given spatial `location`.
"""
function getcell(meta::MetaData, loc)

   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   dx, dy, dz = meta.dx, meta.dy, meta.dz
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells

   @assert xmin < loc[1] ≤ xmax "x coordinate out of bound!"
   @assert ymin < loc[2] ≤ ymax "y coordinate out of bound!"
   @assert zmin < loc[3] ≤ zmax "z coordinate out of bound!"

   # Get cell indices
   indices = floor.(Int, [(loc[1] - xmin)/dx, (loc[2] - ymin)/dy, (loc[3] - zmin)/dz])
   # Get the cell id
   cellid = indices[1] + indices[2] * xcells + indices[3] * xcells * ycells + 1

   # Going through AMR levels as needed
   ilevel = 0
   ncells_lowerlevel = 0

   cellids = meta.cellid
   maxlevel = getmaxamr(meta)

   while ilevel < maxlevel + 1
      if cellid in cellids
         break
      else
         ncells_lowerlevel += 2^(3*ilevel)*(xcells*ycells*zcells)           
         ilevel += 1
         dx *= 0.5; dy *= 0.5; dz *= 0.5

         indices = floor.(Int, [(loc[1] - xmin)/dx, (loc[2] - ymin)/dy, (loc[3] - zmin)/dz])

         cellid = ncells_lowerlevel + indices[1] +
            2^(ilevel)*xcells*indices[2] +
            4^(ilevel)*xcells*ycells*indices[3] + 1
      end
   end
   
   if ilevel == maxlevel + 1
      throw(DomainError(cellid, "CellID does not exist in any AMR level!"))
   end

   cellid
end

"""
    getlevel(meta, cellid) -> Int

Return the AMR level of a given cell ID. Note that this function does not check if the VLSV
file of `meta` actually contains `cellid`: it may be shadowed by refined children.
"""
function getlevel(meta::MetaData, cellid::Integer)
   ncells = meta.xcells*meta.ycells*meta.zcells
   ilevel = 0
   while cellid > 0
      cellid -= 2^(3*ilevel)*ncells
      ilevel += 1
   end
   ilevel - 1 
end

"""
    getmaxamr(meta) -> Int

Find the highest refinement level of a given vlsv file.
"""
function getmaxamr(meta::MetaData)
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells
   maxamr = 0
   cellID = ncells
   while cellID < meta.cellid[end]
      maxamr += 1
      cellID += ncells*8^maxamr
   end

   maxamr
end

"""
    getparent(meta, cellid) -> Int

Return the parent cell ID of given child `cellid`.
"""
function getparent(meta::MetaData, cellid::Integer)
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells

   mylvl = getlevel(meta, cellid)
   parentlvl = mylvl - 1

   if parentlvl < 0
      throw(ArgumentError("$cellid has no parent cell!"))
   else
      # get the first cellid on my level
      cid1st = get1stcell(mylvl, ncells) + 1
      # get my row and column sequence on my level (starting with 0)
      nx = xcells*2^mylvl
      ny = ycells*2^mylvl

      myseq = cellid - cid1st
      ix = myseq % nx
      iz = myseq ÷ (nx*ny)
      iy = (myseq - iz*nx*ny) ÷ nx
      # indexes on the parent level
      ixparent = ix ÷ 2
      iyparent = iy ÷ 2
      izparent = iz ÷ 2
     
      # get the first cellid on parent level
      cid1st -= ncells*8^parentlvl
      # get parent cellid (may not exist!!!)
      parentid = cid1st + izparent*nx*ny÷4 + iyparent*nx÷2 + ixparent
   end
   parentid
end

"""
    getchildren(meta, cellid) -> Vector{Int}

Return direct children of `cellid`.
"""
function getchildren(meta::MetaData, cellid::Integer)
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells

   mylvl = getlevel(meta, cellid)

   # get the first cell ID on the my level
   cid1st = 1
   for i = 0:mylvl-1
      cid1st += ncells*8^i
   end
   # get my row and column sequence on my level (starting with 0)
   nx = xcells*2^mylvl
   ny = ycells*2^mylvl
   
   myseq = cellid - cid1st
   ix = myseq % nx
   iz = myseq ÷ (nx*ny)
   iy = (myseq - iz*nx*ny) ÷ nx

   # get the children sequences on the finer level
   ix *= 2
   iy *= 2
   iz *= 2

   dim = showdimension(meta)
   cid = Vector{Int}(undef, 2^dim)
   # get the first cell ID on the finer level
   cid1st += ncells*8^mylvl
   ix_, iy_ = [ix, ix+1], [iy, iy+1]
   iz_ = zcells != 1 ? [iz, iz+1] : [iz]
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*nx*ny*4 + i[2]*nx*2 + i[1]
   end
   cid
end

"""
    getsiblings(meta, cellid) -> Vector{Int}

Return sibling cells of a given `cellid`, including itself.
"""
function getsiblings(meta::MetaData, cellid::Integer)
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells

   mylvl = getlevel(meta, cellid)

   mylvl == 0 && @error "$cellid is not a child cell!"

   nx = xcells * 2^mylvl
   ny = ycells * 2^mylvl

   # 1st cellid on my level
   cid1st = get1stcell(mylvl, ncells) + 1

   # xyz sequences on my level (starting with 0)
   myseq = cellid - cid1st
   ix = myseq % nx
   iz = myseq ÷ (nx*ny)
   iy = (myseq - iz*nx*ny) ÷ nx
   
   ix1 = iseven(ix) ? ix + 1 : ix - 1
   iy1 = iseven(iy) ? iy + 1 : iy - 1
   iz1 = iseven(iz) ? iz + 1 : iz - 1
   # reorder
   ix, ix1 = minmax(ix, ix1)
   iy, iy1 = minmax(iy, iy1)
   iz, iz1 = minmax(iz, iz1)

   dim = showdimension(meta)
   cid = Vector{Int}(undef, 2^dim)
   ix_, iy_ = [ix, ix1], [iy, iy1]
   iz_ = zcells != 1 ? [iz, iz1] : [iz]
   for (n,i) in enumerate(Iterators.product(ix_, iy_, iz_))
      @inbounds cid[n] = cid1st + i[3]*nx*ny + i[2]*nx + i[1]
   end
   cid
end

"""
    isparent(meta, cellid) -> Bool

Check if `cellid` is a parent cell.
"""
function isparent(meta::MetaData, cellid::Integer)
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells
   amrmax = getmaxamr(meta)

   ncells_accum = get1stcell(amrmax, ncells)

   cellid ∉ meta.cellid && 0 < cellid ≤ ncells_accum 
end

"""
    getcellcoordinates(meta, cellid) -> Vector{Float}

Return a given cell's coordinates.
"""    
function getcellcoordinates(meta::MetaData, cellid::Integer)

   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   cellid -= 1 # for easy divisions

   ncells = xcells*ycells*zcells
   reflevel = 0
   subtraction = ncells * (2^reflevel)^3

   while cellid ≥ subtraction
      cellid -= subtraction
      reflevel += 1
      subtraction *= 8
      xcells *= 2
      ycells *= 2
      zcells *= 2
   end

   indices = zeros(Int, 3)
   indices[1] = cellid % xcells
   indices[2] = cellid ÷ xcells % ycells
   indices[3] = cellid ÷ (xcells*ycells)

   coords = zeros(3)
   coords[1] = xmin + (indices[1] + 0.5) * (xmax - xmin)/xcells
   coords[2] = ymin + (indices[2] + 0.5) * (ymax - ymin)/ycells
   coords[3] = zmin + (indices[3] + 0.5) * (zmax - zmin)/zcells

   coords
end

function isInsideDomain(meta::MetaData, point)
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax

   if xmin < point[1] ≤ xmax && ymin < point[2] ≤ ymax && zmin < point[3] ≤ zmax
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
function getcellinline(meta::MetaData, point1, point2)

   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells

   if !isInsideDomain(meta, point1)
      throw(DomainError(point1, "point location out of bounds!"))
   elseif !isInsideDomain(meta, point2)
      throw(DomainError(point1, "point location out of bounds!"))
   end

   cell_lengths = [(xmax-xmin)/xcells, (ymax-ymin)/ycells, (zmax-zmin)/zcells]

   distances = [0.0]

   cellids = [getcell(meta, point1)]

   coords = point1

   ϵ = eps(Float32)

   p = point1
   unit_vector = @. (point2 - point1) / $hypot(point2 - point1 + ϵ...)

   while true
      cellid = getcell(meta, p)
      amrlvl = getlevel(meta, cellid)

      # Get the max and min cell boundaries
      min_bounds = getcellcoordinates(meta, cellid) -
         0.5 * cell_lengths * 0.5^amrlvl
      max_bounds = min_bounds + cell_lengths

      # Check which face we hit first
      coef_min = (min_bounds - p) ./ unit_vector
      coef_max = (max_bounds - p) ./ unit_vector

      # Negative coefficients indicates the opposite direction.
      @inbounds for i = 1:3
         if unit_vector[i] == 0.0
            coef_min[i] = Inf
            coef_max[i] = Inf
         end
         if coef_min[i] ≤ 0
            coef_min[i] = Inf
         end
         if coef_max[i] ≤ 0
            coef_max[i] = Inf
         end
      end

      # Find the minimum distance from a boundary times a factor
      d = min(minimum(coef_min), minimum(coef_max)) * 1.00001

      coordnew = p + d .* unit_vector

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
    getslicecell(meta, slicelocation, maxamr;
       xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf) -> idlist, indexlist

Find the cell ids `idlist` which are needed to plot a 2d cut through of a 3d mesh, in a
direction with non infinity range at `slicelocation`, and the `indexlist`, which is a
mapping from original order to the cut plane and can be used to select data onto the plane.
"""
function getslicecell(meta::MetaData, slicelocation, maxamr;
   xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf)

   nx, ny, nz = meta.xcells, meta.ycells, meta.zcells
   cellids = meta.cellid # sorted cell IDs

   if !isinf(xmin) && !isinf(xmax)
      minCoord = xmin; maxCoord = xmax; nsize = nx; idim = 1
   elseif !isinf(ymin) && !isinf(ymax)
      minCoord = ymin; maxCoord = ymax; nsize = ny; idim = 2
   elseif !isinf(zmin) && !isinf(zmax)
      minCoord = zmin; maxCoord = zmax; nsize = nz; idim = 3
   else
      @error "Unspecified slice direction!"
   end

   # Find the cut plane index for each refinement level (0-based)
   sliceratio = slicelocation / (maxCoord - minCoord)
   depths = zeros(maxamr+1)
   for i = 0:maxamr
      sliceoffset = floor(Int, sliceratio*nsize*2^i)
      sliceoffset ≤ nsize*2^i || 
         throw(DomainError(sliceoffset, "slice plane index out of bound!"))
      depths[i+1] = sliceoffset
   end

   # Find the ids
   nlen = 0
   ncell = nx*ny*nz
   nCellUptoCurrentLvl = ncell # the number of cells up to refinement level i
   nCellUptoLowerLvl = 0 # the number of cells up to refinement level i-1

   indexlist = Int[]
   idlist = Int[]

   for i = 0:maxamr
      ids = cellids[nCellUptoLowerLvl .< cellids .≤ nCellUptoCurrentLvl]
      ix, iy, iz = getindexes(i, nx, ny, nCellUptoLowerLvl, ids)

      if idim == 1
         coords = ix
      elseif idim == 2
         coords = iy
      elseif idim == 3
         coords = iz
      end

      # Find the needed elements to create the cut and puts the results
      # in the indexlist and the idlist
      elements = coords .== depths[i+1]
      append!(indexlist, (nlen+1:nlen+length(coords))[elements])
      append!(idlist, ids[elements])

      nlen += length(coords)
      nCellUptoLowerLvl = nCellUptoCurrentLvl
      nCellUptoCurrentLvl += ncell*8^(i+1)
   end

   idlist, indexlist
end

"""
    refineslice(meta, idlist, data, maxamr, normal) -> Array

Generate scalar data on the finest refinement level given cellids `idlist` and variable
`data` on the slice perpendicular to `normal`.
"""
function refineslice(meta::MetaData, idlist, data, maxamr, normal)

   nx, ny, nz = meta.xcells, meta.ycells, meta.zcells

   if normal == :x
      dims = [ny, nz] .* 2^maxamr
   elseif normal == :y
      dims = [nx, nz] .* 2^maxamr
   elseif normal == :z
      dims = [nx, ny] .* 2^maxamr
   end

   dpoints = zeros(dims...)

   # Create the plot grid
   ncell = nx*ny*nz
   nCellUptoCurrentLvl = ncell
   nCellUptoLowerLvl = 0

   for i = 0:maxamr
      ids = idlist[nCellUptoLowerLvl .< idlist .≤ nCellUptoCurrentLvl]
      d = data[nCellUptoLowerLvl .< idlist .≤ nCellUptoCurrentLvl]

      ix, iy, iz = getindexes(i, nx, ny, nCellUptoLowerLvl, ids)

      # Get the correct coordinate values and the widths for the plot
      if normal == :x
         a = iy
         b = iz
      elseif normal == :y
         a = ix
         b = iz
      elseif normal == :z
         a = ix
         b = iy
      end

      # Insert the data values into dpoints
      iRange = 0:2^(maxamr - i)-1
      X = [x for x in iRange, _ in iRange]
      Y = [y for _ in iRange, y in iRange]

      coords = Array{Int64,3}(undef, 2, length(a), 2^(2*(maxamr-i)))
      @inbounds for ic = 1:length(a), ir = 1:2^((maxamr-i)*2)
         coords[1,ic,ir] = a[ic]*2^(maxamr - i) + 1 + X[ir]
         coords[2,ic,ir] = b[ic]*2^(maxamr - i) + 1 + Y[ir]
      end

      @inbounds for ic = 1:length(a)
         dpoints[coords[1,ic,:],coords[2,ic,:]] .= d[ic]
      end

      nCellUptoLowerLvl = nCellUptoCurrentLvl
      nCellUptoCurrentLvl += ncell*8^(i+1)
   end

   dpoints
end

"Compute every cell id's x, y and z indexes on the given refinement level (0-based)."
@inline function getindexes(ilevel, nx, ny, nCellUptoLowerLvl, ids)

   slicesize = nx*ny*4^ilevel

   iz = @. (ids - nCellUptoLowerLvl - 1) ÷ slicesize

   # number of ids up to the coordinate z in the refinement level ilevel
   idUpToZ = @. iz*slicesize + nCellUptoLowerLvl

   iy = @. (ids - idUpToZ - 1) ÷ (nx*2^ilevel)
   ix = @. ids - idUpToZ - iy*nx*2^ilevel - 1

   ix, iy, iz
end

function getindexes(ilvl, nx, ny, nCellUptoLowerLvl, id::Int)
   slicesize = nx*ny*4^ilvl
   iz = (id - nCellUptoLowerLvl - 1) ÷ slicesize
   idUpToZ = iz*slicesize + nCellUptoLowerLvl
   iy = (id - idUpToZ - 1) ÷ (nx*2^ilvl)
   ix = id - idUpToZ - iy*nx*2^ilvl - 1
   ix, iy, iz
end

"""
    getnearestcellwithvdf(meta, id) -> Int

Find the nearest spatial cell with VDF saved of a given cell `id` in the file `meta`.
"""
function getnearestcellwithvdf(meta::MetaData, id)
   cells = readmesh(meta.fid, meta.footer, "SpatialGrid", "CELLSWITHBLOCKS")
   isempty(cells) && @error "No distribution saved in $(meta.name)"
   coords = Matrix{Float32}(undef, 3, length(cells))
   for i = 1:length(cells)
      coords[:,i] = getcellcoordinates(meta, cells[i])
   end
   coords_orig = getcellcoordinates(meta, id)
   d2 = sum((coords .- coords_orig).^2, dims=1)
   cells[argmin(d2)[2]]
end


"Return the first cellid - 1 on my level."
function get1stcell(mylevel, ncells)
   cid1st = 0
   for i = 0:mylevel-1
      cid1st += ncells*8^i
   end
   cid1st
end

fillmesh(meta::MetaData, vars::AbstractString) = fillmesh(meta, [vars])

"""
    fillmesh(meta::MetaData, vars; verbose=false)

Fill the DCCRG mesh with quantity of `vars` on all refinement levels.
# Return Arguments
- `celldata::Vector{Vector{Array}}`: data for each variable on each AMR level.
- `vtkGhostType::Array{UInt8}`: cell status (to be completed!). 
"""
function fillmesh(meta::MetaData, vars; verbose=false)

   cellid, maxamr, fid, footer = meta.cellid, meta.maxamr, meta.fid, meta.footer
   nx, ny, nz = meta.xcells, meta.ycells, meta.zcells

   nvarvg = findall(x->!startswith(x, "fg_"), vars)
   nv = length(vars)
   T = Vector{DataType}(undef, nv)
   vsize = Vector{Int}(undef, nv)
   for i = 1:nv
      T[i], _, _, _, vsize[i] = getObjInfo(fid, footer, vars[i], "VARIABLE", "name")
   end

   celldata = [[zeros(T[iv], vsize[iv], nx*2^i, ny*2^i, nz*2^i) for i = 0:maxamr]
      for iv in 1:nv]

   vtkGhostType = [zeros(UInt8, nx*2^i, ny*2^i, nz*2^i) for i = 0:maxamr]

   if maxamr == 0
      for iv = 1:nv
         celldata[iv][1][:] = readvariable(meta, vars[iv])
      end
      return celldata, vtkGhostType
   end

   # Find the ids
   ncell = nx*ny*nz
   nCellUptoCurrentLvl = ncell # the number of cells up to refinement level i
   nCellUptoLowerLvl = 0 # the number of cells up to refinement level i-1

   for ilvl = 0:maxamr
      verbose && @info "scanning AMR level $ilvl..."
      ids = cellid[nCellUptoLowerLvl .< cellid .≤ nCellUptoCurrentLvl]
      
      # indicate the condition of non-existing cells
      idrefined = setdiff(nCellUptoLowerLvl+1:nCellUptoCurrentLvl, ids)

      for id in idrefined
         ix, iy, iz = getindexes(ilvl, nx, ny, nCellUptoLowerLvl, id)
         vtkGhostType[ilvl+1][ix+1,iy+1,iz+1] = 8
      end
      
      if ilvl != maxamr
         for iv in nvarvg
            verbose && @info "reading variable $(vars[iv])..."
            data = readvariable(meta, vars[iv], ids)

            for ilvlup = ilvl:maxamr
               r = 2^(ilvlup-ilvl) # ratio on refined level
               for (ic, id) in enumerate(ids)
                  ix, iy, iz = getindexes(ilvl, nx, ny, nCellUptoLowerLvl, id)
                  for k = 1:r, j = 1:r, i = 1:r
                     @inbounds celldata[iv][ilvlup+1][:,r*ix+i,r*iy+j,r*iz+k] = data[:,ic]
                  end
               end
            end
         end
      else # max amr level
         for (iv, var) = enumerate(vars)
            verbose && @info "reading variable $var..."
            if startswith(var, "fg_")
               celldata[iv][end][:,:,:,:] = readvariable(meta, var)
            else
               data = readvariable(meta, var, ids)
               for (ic, id) in enumerate(ids)
                  ix, iy, iz = getindexes(maxamr, nx, ny, nCellUptoLowerLvl, id)

                  @inbounds celldata[iv][end][:,ix+1,iy+1,iz+1] = data[:,ic]
               end
            end
         end
      end
      nCellUptoLowerLvl = nCellUptoCurrentLvl
      nCellUptoCurrentLvl += ncell*8^(ilvl+1)
   end

   celldata, vtkGhostType
end

"""
    write_vtk(meta::MetaData; vars=[""], ascii=false, verbose=false)
    write_vtk(filename; kwargs...)

Convert VLSV file to VTK format.
# Optional Arguments
- `vars=[""]`: select which variables to convert.
- `ascii=false`: output stored in ASCII or compressed binary format.
- `vti=false`: generate image files on the highest refinement level only.
- `verbose=false`: display logs during conversion.
"""
function write_vtk(meta::MetaData; vars=[""], ascii=false, vti=false, verbose=false)
   nx, ny, nz = meta.xcells, meta.ycells, meta.zcells
   maxamr = meta.maxamr

   append = ascii ? false : true

   filedata = Vector{String}(undef, maxamr+1)
   for i in 1:maxamr+1
      filedata[i] = meta.name[1:end-5]*"_$i.vti"
   end

   if isempty(vars[1])
      vars = showvariables(meta)
      deleteat!(vars, findfirst(x->x=="CellID", vars))
   end

   data, vtkGhostType = fillmesh(meta, vars; verbose)

   if vti
      save_image(meta, meta.name[1:end-4]*"vti", vars, data, vtkGhostType[end], maxamr,
         nx, ny, nz, append)
   else
      # Generate image file on each refinement level
      for i = 1:length(vtkGhostType)
         ghost = vtkGhostType[i]
         save_image(meta, filedata[i], vars, data, ghost, i-1, nx, ny, nz, append)
      end

      # Generate vthb file
      filemeta = meta.name[1:end-4]*"vthb"
      xvthb = XMLDocument()
      xroot = create_root(xvthb, "VTKFile")
      set_attribute(xroot, "type", "vtkOverlappingAMR")
      set_attribute(xroot, "version", "1.1")
      set_attribute(xroot, "byte_order", "LittleEndian") # always the case on x86
      set_attribute(xroot, "header_type", "UInt64")
      xamr = new_child(xroot, "vtkOverlappingAMR")
      origin = @sprintf "%f %f %f" meta.xmin meta.ymin meta.zmin
      set_attribute(xamr, "origin", origin)
      set_attribute(xamr, "grid_description", "XYZ")

      for i = 0:maxamr
         xBlock = new_child(xamr, "Block")
         set_attribute(xBlock, "level", string(i))
         spacing_str = @sprintf "%f %f %f" meta.dx/2^i meta.dy/2^i meta.dz/2^i
         set_attribute(xBlock, "spacing", spacing_str)
         xDataSet = new_child(xBlock, "DataSet")
         set_attribute(xDataSet, "index", "0")
         amr_box = [0, nx*2^i-1, 0, ny*2^i-1, 0, nz*2^i-1]
         box_str = @sprintf "%d %d %d %d %d %d" amr_box...
         set_attribute(xDataSet, "amr_box", box_str)
         set_attribute(xDataSet, "file", filedata[i+1])
      end

      save_file(xvthb, filemeta)
   end

   return
end

write_vtk(filename; kwargs...) = write_vtk(readmeta(filename); kwargs...)

"""
    save_image(meta::MetaData, file, vars, data, vtkGhostType, level, nx, ny, nz)

Save `data` of name `vars` at AMR `level` into VTK image file of name `file`.
# Arguments
`file::String`: output file name.
`vars::Vector{String}`: variable names to be saved.
`data::Vector{Vector}`: data for all the variables on each refinement level.
`vtkGhostType::Array{UInt8}`: array for visibility control.
`level::Int`: refinement level (0-based).
`nx, ny, nz`: original mesh sizes.
`ascii=false`: save output in ASCII or binary format.
`append=true`: determines whether to append data at the end of file or do in-block writing.
"""
function save_image(meta::MetaData, file, vars, data, vtkGhostType, level, nx, ny, nz,
   ascii=false, append=true)
   origin = (meta.xmin, meta.ymin, meta.zmin)
   ratio = 2^level
   spacing = (meta.dx / ratio, meta.dy / ratio, meta.dz / ratio)

   vtk = vtk_grid(file, nx*ratio+1, ny*ratio+1, nz*ratio+1; origin, spacing, append, ascii)

   for (iv, var) in enumerate(vars)
      vtk[var, VTKCellData()] = data[iv][level+1]
   end

   vtk["vtkGhostType", VTKCellData()] = vtkGhostType

   vtk_save(vtk)
end


"""
    compare(filename1, filename2, tol=1e-4) -> Bool

Check if two VLSV files are approximately identical.
"""
function compare(f1, f2, tol::AbstractFloat=1e-4)
   # 1st sanity check: minimal filesize difference
   if abs(filesize(f1) - filesize(f2)) / filesize(f2) > 1e-2
      return false
   end

   meta1 = readmeta(f1)
   meta2 = readmeta(f2)
   varnames = showvariables(meta1)
   strskip = r"CellID|rank|blocks"
   deleteat!(varnames, findall(x->endswith(x, strskip), varnames))

   isIdentical = true
   for vname in varnames
      v1 = readvariable(meta1, vname)
      v2 = readvariable(meta2, vname)

      s1, s2 = sum(v1), sum(v2)
      if abs(s1 - s2) > tol * abs(s1) && abs(s1 - s2) > tol * abs(s2)
         isIdentical = false
         println("$vname is quite different!")
         break
      end
   end
   close(meta1.fid)
   close(meta2.fid)
   return isIdentical
end