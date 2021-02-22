# Utility functions for processing VLSV data.

import LinearAlgebra: norm, dot

const qₑ = -1.60217662e-19  # electron charge, [C]
const mₑ = 9.10938356e-31   # electron mass, [kg]
const qᵢ = 1.60217662e-19   # proton mass, [C]
const mᵢ = 1.673557546e-27  # proton mass, [kg]
const c  = 3e8              # speed of light, [m/s]
const μ₀ = 4π*1e-7          # Vacuum permeability, [H/m]
const Re = 6.371e6          # Earth radius, [m]

export get_cellid, getSliceCellID, get_amr_level, get_max_amr_level,
   refine_data, get_cell_coordinates, get_cell_in_line,
   getNearestCellWithVspace, compare

"""
    get_cell_id(meta::MetaData, location)

Return cell ID containing the given spatial location.
"""
function get_cellid(meta, loc)

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
   maxlevel = get_max_amr_level(meta)

   while ilevel < maxlevel + 1
      if cellid in cellids
         break
      else
         ncells_lowerlevel += 2^(3*ilevel)*(xcells*ycells*zcells)           
         ilevel += 1
         dx /= 2; dy /= 2; dz /= 2

         indices = floor.(Int, [(loc[1] - xmin)/dx, (loc[2] - ymin)/dy, (loc[3] - zmin)/dz])

         cellid = ncells_lowerlevel + indices[1] +
            2^(ilevel)*xcells*indices[2] +
            4^(ilevel)*xcells*ycells*indices[3] + 1
      end
   end
   
   if ilevel == maxlevel + 1
      @error "CellID does not exist in any AMR level"
   end

   return cellid
end

"""
    get_amr_level(meta::MetaData, cellid)

Return the AMR level of a given cell id.
"""
function get_amr_level(meta, cellid)

   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells

   ilevel = 0
   while cellid > 0
      cellid -= 2^(3*(ilevel))*(xcells*ycells*zcells)
      ilevel += 1
   end
   return ilevel - 1 
end

"""
    get_max_amr_level(meta::MetaData)

Find the highest refinement level.
"""
function get_max_amr_level(meta)
   xcells, ycells, zcells= meta.xcells, meta.ycells, meta.zcells
   ncells = xcells*ycells*zcells
   maxreflevel = 0
   cellID = ncells
   while cellID < meta.cellid[end]
      maxreflevel += 1
      cellID += Int(ncells*8^maxreflevel)
   end

   maxreflevel
end

"""
    get_cell_coordinates(meta::MetaData, cellid)

Return a given cellid's coordinates.
"""    
function get_cell_coordinates(meta, cellid)

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

   return coords
end

function isInsideDomain(meta, point)
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax

   if xmin < point[1] ≤ xmax && ymin < point[2] ≤ ymax && zmin < point[3] ≤ zmax
      return true
   else
      return false
   end
end

"""
    get_cell_in_line(meta::MetaData, point1, point2)

Returns cell IDs, distances and coordinates for every cell in a line between two
given points. May be improved later with preallocation!
"""
function get_cell_in_line(meta, point1, point2)

   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   xcells, ycells, zcells = meta.xcells, meta.ycells, meta.zcells
   dx, dy, dz = meta.dx, meta.dy, meta.dz

   if !isInsideDomain(meta, point1)
      @error "point1 in get_cell_in_line out of bounds!"
   elseif !isInsideDomain(meta, point2)
      @error "point2 in get_cell_in_line out of bounds!"
   end

   cell_lengths = [(xmax-xmin)/xcells, (ymax-ymin)/ycells, (zmax-zmin)/zcells]

   distances = [0.0]

   cellids = [get_cellid(meta, point1)]

   coords = point1

   ϵ = eps(Float32)

   p = point1
   unit_vector = @. (point2 - point1) / $norm(point2 - point1 + ϵ)

   while true
      cellid = get_cellid(meta, p)
      amrlvl = get_amr_level(meta, cellid)

      # Get the max and min cell boundaries
      min_bounds = get_cell_coordinates(meta, cellid) -
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
      d = min(minimum(coef_min), minimum(coef_max)) * 1.0001

      coordnew = p + d .* unit_vector
      cellidnew = get_cellid(meta, coordnew)

      push!(cellids, cellidnew)
      coords = hcat(coords, coordnew)
      push!(distances, norm(coordnew - point1))

      p = coordnew

      dot(point2-p, unit_vector) ≥ 0 || break
   end

   return cellids, distances, coords
end

"""
    getSliceCellID(meta, slicelocation, maxreflevel; xmin=-Inf, xmax=Inf,
       ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf)

Find the cell ids which are needed to plot a 2d cut through of a 3d mesh.
"""
function getSliceCellID(meta, slicelocation, maxreflevel;
   xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, zmin=-Inf, zmax=Inf)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells
   cellids = meta.cellid # sorted cell IDs

   if !isinf(xmin) && !isinf(xmax)
      minCoord = xmin; maxCoord = xmax; nsize = xsize; idim = 1
   elseif !isinf(ymin) && !isinf(ymax)
      minCoord = ymin; maxCoord = ymax; nsize = ysize; idim = 2
   elseif !isinf(zmin) && !isinf(zmax)
      minCoord = zmin; maxCoord = zmax; nsize = zsize; idim = 3
   else
      @error "Unspecified slice direction!"
   end

   # Find the cut plane index for each refinement level
   sliceratio = slicelocation/(maxCoord-minCoord)
   depths = zeros(maxreflevel+1)
   for i = 0:maxreflevel
      sliceoffset = floor(Int32, sliceratio*nsize*2^i) + 1
      sliceoffset ≤ nsize*2^i || @error "slice plane index out of bound!"
      depths[i+1] = sliceoffset
   end

   # Find the ids
   nlen = 0
   ncell = xsize*ysize*zsize
   nCellUptoCurrentLvl = ncell # the number of cells up to refinement level i
   nCellUptoLowerLvl = 0 # the number of cells up to refinement level i-1

   indexlist = Int[]
   idlist = Int[]

   for i = 0:maxreflevel
      ids = cellids[nCellUptoLowerLvl .< cellids .≤ nCellUptoCurrentLvl]

      # Compute every cell ids' x, y and z indexes
      z = @. (ids - nCellUptoLowerLvl - 1) ÷ (xsize*ysize*4^i) + 1

      # number of ids up to the coordinate z in the refinement level i
      idUpToZ = @. (z-1)*xsize*ysize*4^i + nCellUptoLowerLvl

      y = @. (ids - idUpToZ - 1) ÷ (xsize*2^i) + 1
      x = @. ids - idUpToZ - (y-1)*xsize*2^i

      if idim == 1
         coords = x
      elseif idim == 2
         coords = y
      elseif idim == 3
         coords = z
      end

      # Find the needed elements to create the cut and puts the results in the
      # indexlist and the idlist
      elements = coords .== depths[i+1]
      append!(indexlist, (nlen+1:nlen+length(coords))[elements])
      append!(idlist, ids[elements])

      nlen += length(coords)
      nCellUptoLowerLvl = nCellUptoCurrentLvl
      nCellUptoCurrentLvl += ncell*8^(i+1)
   end

   return idlist, indexlist
end

"""
    refine_data(meta, idlist, data, maxreflevel, normal)

Generate scalar data on the finest refinement level.
"""
function refine_data(meta, idlist, data, maxreflevel, normal)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells

   if normal == :x
      dims = [ysize, zsize] .* 2^maxreflevel
   elseif normal == :y
      dims = [xsize, zsize] .* 2^maxreflevel
   elseif normal == :z
      dims = [xsize, ysize] .* 2^maxreflevel
   end

   dpoints = zeros(dims...)

   # Create the plot grid
   ncell = xsize*ysize*zsize
   nCellUptoCurrentLvl = ncell
   nCellUptoLowerLvl = 0

   for i = 0:maxreflevel
      ids = idlist[nCellUptoLowerLvl .< idlist .≤ nCellUptoCurrentLvl]

      d = data[nCellUptoLowerLvl .< idlist .≤ nCellUptoCurrentLvl]

      # Compute every cell ids' x, y and z indexes on this refinement level
      z = @. (ids - nCellUptoLowerLvl - 1) ÷ (xsize*ysize*4^i) + 1

      # number of ids up to the coordinate z in the refinement level i
      idUpToZ = @. (z-1)*xsize*ysize*4^i + nCellUptoLowerLvl

      y = @. (ids - idUpToZ - 1) ÷ (xsize*2^i) + 1
      x = @. ids - idUpToZ - (y-1)*xsize*2^i

      # Get the correct coordinate values and the widths for the plot
      if normal == :x
         a = y
         b = z
      elseif normal == :y
         a = x
         b = z
      elseif normal == :z
         a = x
         b = y
      end

      # Insert the data values into dpoints
      iRange = 0:2^(maxreflevel - i)-1
      X = [x for x in iRange, _ in iRange]
      Y = [y for _ in iRange, y in iRange]

      coords = Array{Int64,3}(undef, 2, length(a), 2^(2*(maxreflevel-i)))
      @inbounds for ic = 1:length(a), ir = 1:2^((maxreflevel-i)*2)
         coords[1,ic,ir] = (a[ic] - 1)*2^(maxreflevel - i) + 1 + X[ir]
         coords[2,ic,ir] = (b[ic] - 1)*2^(maxreflevel - i) + 1 + Y[ir]
      end

      @inbounds for ic = 1:length(a)
         dpoints[coords[1,ic,:],coords[2,ic,:]] .= d[ic]
      end

      nCellUptoLowerLvl = nCellUptoCurrentLvl
      nCellUptoCurrentLvl += ncell*8^(i+1)
   end

   return dpoints
end

"Find the nearest spatial cell with f saved of a given cell `id`."
function getNearestCellWithVspace(meta, id)
   cells = Vlasiator.read_mesh(meta.fid, meta.footer, "SpatialGrid",
      "CELLSWITHBLOCKS")
   coords = Matrix{Float32}(undef, 3, length(cells))
   for i = 1:length(cells)
      coords[:,i] = get_cell_coordinates(meta, cells[i])
   end
   coords_orig = get_cell_coordinates(meta, id)
   d2 = sum((coords .- coords_orig).^2, dims=1)
   cells[argmin(d2)[2]]
end

"Compare if two VLSV files are identical."
function compare(f1, f2)
   meta1 = read_meta(f1)
   meta2 = read_meta(f2)
   varnames = show_variables(meta1)
   list_nocheck = ("CellID", "vg_rank")
   deleteat!(varnames, findall(x->x in list_nocheck, varnames))
   isIdentical = true
   for vname in varnames
      v1 = read_variable(meta1, vname)
      v2 = read_variable(meta2, vname)
      v1 != v2 && (isIdentical = false; break)
   end
   close(meta1.fid)
   close(meta2.fid)
   return isIdentical
end