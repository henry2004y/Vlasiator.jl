# VLSV reader in Julia

include("vlsvvariables.jl")

using LightXML

export MetaData, VarInfo
export hasvariable, hasparameter, hasname, hasvdf
export readmeta, readvariable, readparameter, readvariablemeta, readvcells
export showvariables, showdimension, getvcellcoordinates

"Mesh size information."
struct MeshInfo
   vxblocks::Int64
   vyblocks::Int64
   vzblocks::Int64
   vxblock_size::Int64
   vyblock_size::Int64
   vzblock_size::Int64
   vxmin::Float64
   vymin::Float64
   vzmin::Float64
   vxmax::Float64
   vymax::Float64
   vzmax::Float64
   dvx::Float64
   dvy::Float64
   dvz::Float64
end

"""
Variable metadata from the vlsv footer.
"""
struct VarInfo
   "unit of the variable as a string"
   unit::String
   "unit of the variable as a LaTeX-formatted string"
   unitLaTeX::LaTeXString
   "description of the variable as a LaTeX-formatted string"
   variableLaTeX::LaTeXString
   "conversion factor to SI units as a string"
   unitConversion::String
end

"Meta data declaration."
struct MetaData
   name::AbstractString
   fid::IOStream
   footer::XMLElement
   cellid::Vector{UInt64}  # sorted cell IDs
   cellIndex::Vector{Int64}
   maxamr::Int64
   xcells::Int64
   ycells::Int64
   zcells::Int64
   xblock_size::Int64
   yblock_size::Int64
   zblock_size::Int64
   xmin::Float64
   ymin::Float64
   zmin::Float64
   xmax::Float64
   ymax::Float64
   zmax::Float64
   dx::Float64
   dy::Float64
   dz::Float64
   meshes::Dict{String, MeshInfo}
   populations::Vector{String}
end


function Base.show(io::IO, meta::MetaData)
   println(io, "filename = ", meta.name)
   dim = showdimension(meta)
   println(io, "dimension: $dim")
   println(io, "maximum AMR level: $(meta.maxamr)")
   println(io, "contains VDF: $(hasvdf(meta))")
   vars = showvariables(meta)
   print(io, "variables: ")
   println(io, vars)
end

function Base.show(io::IO, s::VarInfo)
   println(io, "var in LaTeX: ", s.variableLaTeX)
   println(io, "unit: ", s.unit)
   println(io, "unit in LaTeX: ", s.unitLaTeX)
   println(io, "unit conversion = ", s.unitConversion)
end

"Return the xml footer of vlsv."
function getfooter(fid)

   # First 8 bytes indicate big-endian or else
   endianness_offset = 8
   seek(fid, endianness_offset)
   # Obtain the offset of the XML file
   uint64_byte_amount = 8
   offset = read(fid, UInt64)
   seek(fid, offset)
   xmldata = read(fid, String)
   xmldoc  = parse_string(xmldata)
   footer = root(xmldoc)

end


"Return size and type information for the object."
function getObjInfo(fid, footer, name, tag, attr)

   arraysize = 0
   datasize = 0
   datatype = ""
   vectorsize = 0
   variable_offset = 0

   isFound = false
   
   for varinfo in footer[tag]
      if attribute(varinfo, attr) == name
         arraysize = parse(Int, attribute(varinfo, "arraysize"))
         datasize = parse(Int, attribute(varinfo, "datasize"))
         datatype = attribute(varinfo, "datatype")
         vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
         variable_offset = parse(Int, content(varinfo))
         isFound = true
         break
      end
   end

   if !isFound @error "unknown variable $(name)!" end

   seek(fid, variable_offset)

   if datatype == "float" && datasize == 4
      T = Float32
   elseif datatype == "float" && datasize == 8
      T = Float64
   elseif datatype == "int" && datasize == 4
      T = Int32
   elseif datatype == "int" && datasize == 8
      T = Int64
   elseif datatype == "uint" && datasize == 4
      T = UInt32
   elseif datatype == "uint" && datasize == 8
      T = UInt64
   end

   T::DataType, variable_offset, arraysize, datasize, vectorsize
end

"Return vector data from vlsv file."
function readvector(fid, footer, name, tag)

   T, _, arraysize, _, vectorsize = getObjInfo(fid, footer, name, tag, "name")

   if vectorsize == 1
      w = Vector{T}(undef, arraysize)
   else
      w = Array{T,2}(undef, vectorsize, arraysize)
   end

   read!(fid, w)

   w
end


"""
    readmeta(filename; verbose=false) -> MetaData

Return MetaData from a vlsv file.
"""
function readmeta(filename::AbstractString; verbose=false)

   fid = open(filename, "r")

   footer = getfooter(fid)

   meshName = "SpatialGrid"

   cellid = readvector(fid, footer, "CellID", "VARIABLE")

   cellIndex = sortperm(cellid)

   bbox = readmesh(fid, footer, meshName, "MESH_BBOX")

   nodeCoordsX = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_X")
   nodeCoordsY = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_Y")
   nodeCoordsZ = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_Z")
  
   xcells, ycells, zcells = bbox[1:3]
   xblock_size, yblock_size, zblock_size = bbox[4:6]
   xmin, ymin, zmin = nodeCoordsX[1], nodeCoordsY[1], nodeCoordsZ[1]
   xmax, ymax, zmax = nodeCoordsX[end], nodeCoordsY[end], nodeCoordsZ[end]

   dx = (xmax - xmin) / xcells
   dy = (ymax - ymin) / ycells
   dz = (zmax - zmin) / zcells

   meshes = Dict{String,MeshInfo}()

   # Find all populations by the BLOCKIDS tag
   populations = String[]

   for varinfo in footer["BLOCKIDS"]

      if has_attribute(varinfo, "name")
         # New style vlsv file with bounding box
         popname = attribute(varinfo, "name")

         bbox = readmesh(fid, footer, popname, "MESH_BBOX")

         nodeCoordsX = readmesh(fid, footer, popname, "MESH_NODE_CRDS_X")   
         nodeCoordsY = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Y")   
         nodeCoordsZ = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Z")   
         vxblocks, vyblocks, vzblocks = bbox[1:3]
         vxblock_size, vyblock_size, vzblock_size = bbox[4:6]
         vxmin = nodeCoordsX[1]
         vymin = nodeCoordsY[1]
         vzmin = nodeCoordsZ[1]
         vxmax = nodeCoordsX[end]
         vymax = nodeCoordsY[end]
         vzmax = nodeCoordsZ[end]
         dvx = (vxmax - vxmin) / vxblocks / vxblock_size
         dvy = (vymax - vymin) / vyblocks / vyblock_size
         dvz = (vzmax - vzmin) / vzblocks / vzblock_size
      else
         popname = "avgs"

         if "vxblocks_ini" in attribute.(footer["PARAMETER"],"name") 
            # Old vlsv files where the mesh is defined with parameters
            vxblocks = readparameter(fid, footer, "vxblocks_ini")
            vyblocks = readparameter(fid, footer, "vyblocks_ini")
            vzblocks = readparameter(fid, footer, "vzblocks_ini")
            vxblock_size = 4
            vyblock_size = 4
            vzblock_size = 4
            vxmin = readparameter(fid, footer, "vxmin")
            vymin = readparameter(fid, footer, "vymin")
            vzmin = readparameter(fid, footer, "vzmin")
            vxmax = readparameter(fid, footer, "vxmax")
            vymax = readparameter(fid, footer, "vymax")
            vzmax = readparameter(fid, footer, "vzmax")
            dvx = (vxmax - vxmin) / vxblocks / vxblock_size
            dvy = (vymax - vymin) / vyblocks / vyblock_size
            dvz = (vzmax - vzmin) / vzblocks / vzblock_size
         else
            # No velocity space info, e.g., file not written by Vlasiator 
            vxblocks, vyblocks, vzblocks = 0, 0, 0
            vxblock_size, vyblock_size, vzblock_size = 4, 4, 4
            vxmin, vymin, vzmin = 0.0, 0.0, 0.0
            vxmax, vymax, vzmax = 0.0, 0.0, 0.0
            dvx, dvy, dvz = 1.0, 1.0, 1.0
         end
      end

      # Update list of active populations
      if popname ∉ populations 
         push!(populations, popname)
      end

      # Create a new MeshInfo object for this population
      popMesh = MeshInfo(vxblocks, vyblocks, vzblocks, 
         vxblock_size, vyblock_size, vzblock_size,
         vxmin, vymin, vzmin, vxmax, vymax, vzmax,
         dvx, dvy, dvz)

      meshes[popname] = popMesh

      if verbose
         @info "Found population " * popname
      end
   end
   
   # Obtain maximum refinement level
   ncells = xcells*ycells*zcells
   maxamr = 0
   cid = ncells
   while cid < cellid[cellIndex[end]]
      maxamr += 1
      cid += ncells*8^maxamr
   end

   #close(fid) # Is it safe not to close it?

   meta = MetaData(filename, fid, footer, cellid[cellIndex], cellIndex, maxamr,
      xcells, ycells, zcells, xblock_size, yblock_size, zblock_size,
      xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, meshes, populations)
end


"""
    readvariablemeta(meta, var) -> VarInfo

Return VarInfo about `var` in the vlsv file linked to `meta`.           
"""
function readvariablemeta(meta, var)

   unit = ""
   unitLaTeX = ""
   variableLaTeX = ""
   unitConversion = ""

   # Force lowercase
   var = lowercase(var)

   # Get population and variable names from data array name 
   if occursin("/", var)
      popname, varname = split(var, "/")
   else
      popname = "pop"
      varname = var
   end

   if hasvariable(meta, var)
      if varname[1:3] == "vg_" || varname[1:3] == "fg_"
         # For Vlasiator 5 vlsv files, metadata is included

         for varinfo in meta.footer["VARIABLE"]
            if attribute(varinfo, "name") == var
               unit = attribute(varinfo, "unit")
               unitLaTeX = attribute(varinfo, "unitLaTeX")
               variableLaTeX = attribute(varinfo, "variableLaTeX")
               unitConversion = attribute(varinfo, "unitConversion") 
            end
         end

      elseif var in keys(units_predefined)
         unit = units_predefined[var]
         variableLaTeX = latex_predefined[var]
         unitLaTeX = latexunits_predefined[var]
      end
   end

   VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

"Return mesh related variable."
function readmesh(fid, footer, typeMesh, varMesh)

   T, variable_offset, arraysize, datasize, vectorsize = 
      getObjInfo(fid, footer, typeMesh, varMesh, "mesh")

   w = Vector{T}(undef, arraysize)
   read!(fid, w)

   w
end

"""
    readvariable(meta, var, sorted=true) -> Array

Return variable value from the vlsv file. By default `sorted=true`, which means that for
DCCRG grid the variables are sorted by cell ID.
"""
function readvariable(meta::MetaData, var::AbstractString, sorted::Bool=true)

   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
      return data
   end

   data = readvector(meta.fid, meta.footer, var, "VARIABLE")

   if startswith(var, "fg_") # fsgrid
      bbox = readmesh(meta.fid, meta.footer, "fsgrid", "MESH_BBOX")
      # Determine fsgrid domain decomposition
      nIORanks = readparameter(meta, "numWritingRanks")

      if ndims(data) > 1
         orderedData = zeros(Float32, size(data,1), bbox[1:3]...)
      else
         orderedData = zeros(Float32, bbox[1:3]...)
      end

      currentOffset = 1
      fgDecomposition = getDomainDecomposition(bbox[1:3], nIORanks)

      @inbounds for i = 0:nIORanks-1
         x = i ÷ fgDecomposition[3] ÷ fgDecomposition[2]
         y = i ÷ fgDecomposition[3] % fgDecomposition[2]
         z = i % fgDecomposition[3]

         lsize = calcLocalSize.(bbox[1:3], fgDecomposition, [x,y,z])
         lstart = calcLocalStart.(bbox[1:3], fgDecomposition, [x,y,z])
         lend = @. lstart + lsize - 1

         totalSize = prod(lsize)

         # Reorder data
         if ndims(data) > 1
            ldata = data[:,currentOffset:currentOffset+totalSize-1]
            ldata = reshape(ldata, size(data,1), lsize...)

            orderedData[:,lstart[1]:lend[1],lstart[2]:lend[2],lstart[3]:lend[3]] = ldata
         else
            ldata = data[currentOffset:currentOffset+totalSize-1]
            ldata = reshape(ldata, lsize...)

            orderedData[lstart[1]:lend[1],lstart[2]:lend[2],lstart[3]:lend[3]] = ldata
         end

         currentOffset += totalSize
      end
      data = dropdims(orderedData, dims=(findall(size(orderedData) .== 1)...,))
   elseif sorted # dccrg grid
      if ndims(data) == 1
         data = data[meta.cellIndex]
      elseif ndims(data) == 2
         data = data[:, meta.cellIndex]
      end
   end
   return data
end

"""
    readvariable(meta, var, ids) -> Array

Read a variable `var` in a collection of cells `ids`.
"""
function readvariable(meta::MetaData, var::AbstractString, ids)

   @assert !startswith(var, "fg_") "Currently does not support reading fsgrid!"

   if isempty(ids)
      w = readvariable(meta, var, false)
      return [w]
   end

   T, variable_offset, arraysize, datasize, vectorsize = 
      getObjInfo(meta.fid, meta.footer, var, "VARIABLE", "name")

   cellids = readvariable(meta, "CellID", false)
   rOffsets = [(findfirst(x->x==i, cellids)-1)*datasize*vectorsize for i in ids]

   v = Array{T}(undef, vectorsize, length(ids))

   for (i, r) in enumerate(rOffsets)
      loc = variable_offset + r
      seek(meta.fid, loc)
   
      w = Vector{T}(undef, vectorsize)
      read!(meta.fid, w)
      v[:,i] = w
   end

   return v
end

# Optimize decomposition of this grid over the given number of processors.
# Reference: fsgrid.hpp
function getDomainDecomposition(globalsize, nprocs)
   domainDecomp = [1, 1, 1]
   procBox = [0.0, 0.0, 0.0]
   minValue = Inf

   for i = 1:min(nprocs, globalsize[1])
      procBox[1] = max(globalsize[1]/i, 1)

      for j = 1:min(nprocs, globalsize[2])
         i * j > nprocs && break

         procBox[2] = max(globalsize[2]/j, 1)

         for k = 1:min(nprocs, globalsize[2])
            i * j * k > nprocs && continue

            procBox[3] = max(globalsize[3]/k, 1)
         
            nyz = i > 1 ? procBox[2] * procBox[3] : 0
            nzx = j > 1 ? procBox[1] * procBox[3] : 0
            nxy = k > 1 ? procBox[1] * procBox[2] : 0

            v = 10*procBox[1]*procBox[2]*procBox[3] + nyz + nzx + nxy            

            if i * j * k == nprocs && v < minValue
               minValue = v
               domainDecomp[1:3] = [i, j, k]
            end
         end
      end
   end

   domainDecomp
end

function calcLocalStart(globalCells, nprocs, lcells)
   ncells = globalCells ÷ nprocs
   remainder = globalCells % nprocs
   lstart = lcells < remainder ? 
      lcells*(ncells+1) + 1 :
      lcells*ncells + remainder + 1
end

function calcLocalSize(globalCells, nprocs, lcells)
   ncells = globalCells ÷ nprocs
   remainder = globalCells % nprocs
   lsize = lcells < remainder ? ncells + 1 : ncells
end

"""
    hasvariable(meta, var) -> Bool

Check if the VLSV file contains a variable.
"""
hasvariable(meta::MetaData, var) = hasname(meta.footer, "VARIABLE", var)

"""
    readparameter(meta, param)

Return the parameter value from vlsv file.
"""
readparameter(meta::MetaData, param) = readparameter(meta.fid, meta.footer, param)

function readparameter(fid, footer, param)
   
   T, _, _, _, _ = getObjInfo(fid, footer, param, "PARAMETER", "name")

   p = read(fid, T)
end

"""
    hasparameter(meta, param) -> Bool

Check if vlsv file contains a parameter.
"""
hasparameter(meta::MetaData, param) = hasname(meta.footer, "PARAMETER", param)

"Check if the XMLElement `elem` contains a `tag` with `name`."
function hasname(elem, tag, name)
   isFound = false
   
   for varinfo in elem[tag]
      attribute(varinfo, "name") == name && (isFound = true)
      isFound && break
   end
   
   isFound
end

"""
    showvariables(meta) -> Vector{String}

Return all variable names in the VLSV file.
"""
function showvariables(meta::MetaData)
   nVar = length(meta.footer["VARIABLE"])
   vars = Vector{String}(undef, nVar)
   for i in 1:nVar
      vars[i] = attribute(meta.footer["VARIABLE"][i], "name")
   end
   vars
end

"""
    showdimension(meta) -> Int

Return the dimension of data.
"""
showdimension(meta::MetaData) = count(>(1), [meta.xcells, meta.ycells, meta.zcells])

"""
    hasvdf(meta) -> Bool

Check if VLSV file contains VDF.
"""
function hasvdf(meta::MetaData) 
   cells = readmesh(meta.fid, meta.footer, "SpatialGrid", "CELLSWITHBLOCKS")
   if isempty(cells)
      return false
   else
      return true
   end
end

"""
    readvcells(meta, cellid; pop="proton")

Read velocity cells from a spatial cell of ID `cellid`, and return a map of velocity cell
ids and corresponding value.
"""
function readvcells(meta, cellid; pop="proton")

   fid, footer, vmesh = meta.fid, meta.footer, meta.meshes[pop]
   nblockx = vmesh.vxblock_size
   nblocky = vmesh.vyblock_size
   nblockz = vmesh.vzblock_size
   bsize = nblockx * nblocky * nblockz

   cellsWithVDF = readvector(fid, footer, pop, "CELLSWITHBLOCKS")
   nblock_C = readvector(fid, footer, pop, "BLOCKSPERCELL")

   nblock_C_offsets = zeros(Int, length(cellsWithVDF))
   nblock_C_offsets[2:end] = cumsum(nblock_C[1:end-1])

   # Check if cells have vspace stored
   if cellid ∈ cellsWithVDF
      cellWithVDFIndex = findfirst(x->x==cellid, cellsWithVDF)
   else
      @error "The input cell does not have velocity distribution!"
   end

   # Navigate to the correct position
   offset = nblock_C_offsets[cellWithVDFIndex]
   nblocks = nblock_C[cellWithVDFIndex]

   # Read in avgs
   varinfo = footer["BLOCKVARIABLE"][1]
   arraysize = parse(Int, attribute(varinfo, "arraysize"))
   datasize = parse(Int, attribute(varinfo, "datasize"))
   datatype = attribute(varinfo, "datatype")
   vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
   variable_offset = parse(Int, content(varinfo))

   # Navigate to the correct position
   offset_data = offset * vectorsize * datasize + variable_offset

   seek(fid, offset_data)

   if datatype == "float" && datasize == 4
      Tavg = Float32
   elseif datatype == "float" && datasize == 8
      Tavg = Float64
   end

   data = Array{Tavg,2}(undef, vectorsize, nblocks)
   read!(fid, data)

   # Read in block IDs
   varinfo = footer["BLOCKIDS"][1]
   arraysize = parse(Int, attribute(varinfo, "arraysize"))
   datasize = parse(Int, attribute(varinfo, "datasize"))
   datatype = attribute(varinfo, "datatype")
   vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
   variable_offset = parse(Int, content(varinfo))

   # Navigate to the correct position
   offset_f = offset * vectorsize * datasize + variable_offset

   seek(fid, offset_f)

   if datatype == "uint" && datasize == 4
      T = UInt32
   elseif datatype == "uint" && datasize == 8
      T = UInt64
   end

   blockIDs = Vector{T}(undef, nblocks*vectorsize)
   read!(fid, blockIDs)

   # Velocity cell IDs and corresponding avg distributions
   vcellids = zeros(Int, bsize*nblocks)
   vcellf = zeros(Tavg, bsize*nblocks)

   vcellid_local = [i + nblockx*j + nblockx*nblocky*k
      for i in 0:nblockx-1, j in 0:nblocky-1, k in 0:nblockz-1]

   @inbounds for i in 1:nblocks
      vblockid = blockIDs[i]
      for j = 1:bsize
         vcellids[(i-1)*bsize+j] = vcellid_local[j] + bsize*vblockid
         vcellf[(i-1)*bsize+j] = data[j,i]
      end
   end
   vcellids, vcellf
end


"""
    getvcellcoordinates(meta, vcellids, pop="proton")

Return velocity cells' coordinates of population `pop` and id `vcellids`.
"""
function getvcellcoordinates(meta, vcellids; pop="proton")

   vmesh = meta.meshes[pop]
   bsize = vmesh.vxblock_size * vmesh.vyblock_size * vmesh.vzblock_size
   blockid = @. vcellids ÷ bsize
   # Get block coordinates
   blockIndX = @. blockid % (vmesh.vxblocks)
   blockIndY = @. blockid ÷ vmesh.vxblocks % vmesh.vyblocks
   blockIndZ = @. blockid ÷ (vmesh.vxblocks * vmesh.vyblocks)
   blockCoordX = @. blockIndX * vmesh.dvx * vmesh.vxblock_size + vmesh.vxmin
   blockCoordY = @. blockIndY * vmesh.dvy * vmesh.vyblock_size + vmesh.vymin
   blockCoordZ = @. blockIndZ * vmesh.dvz * vmesh.vzblock_size + vmesh.vzmin
   # Get cell indices
   cellids = @. vcellids % bsize
   cellidx = @. cellids % vmesh.vxblock_size
   cellidy = @. cellids ÷ vmesh.vxblock_size % vmesh.vyblock_size
   cellidz = @. cellids ÷ (vmesh.vxblock_size * vmesh.vyblock_size)
   # Get cell coordinates
   cellCoords = Matrix{Float32}(undef, 3, length(cellids))
   @inbounds for i = 1:length(cellids)
      cellCoords[1,i] = blockCoordX[i] + (cellidx[i] + 0.5) * vmesh.dvx
      cellCoords[2,i] = blockCoordY[i] + (cellidy[i] + 0.5) * vmesh.dvy
      cellCoords[3,i] = blockCoordZ[i] + (cellidz[i] + 0.5) * vmesh.dvz
   end
   cellCoords
end