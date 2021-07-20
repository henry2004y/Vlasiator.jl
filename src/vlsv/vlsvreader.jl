# VLSV reader in Julia

include("vlsvvariables.jl")

using LightXML, FLoops

export MetaData, VarInfo
export hasvariable, hasparameter, hasname, hasvdf
export load, readvariable, readparameter, readvariablemeta, readvcells
export ndims, getvcellcoordinates

"Velocity mesh information."
struct VMeshInfo
   "number of velocity blocks"
   vblocks::Vector{Int64}
   vblock_size::Vector{Int64}
   vmin::Vector{Float64}
   vmax::Vector{Float64}
   dv::Vector{Float64}
end

"Variable metadata from the vlsv footer."
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

"VLSV meta data."
struct MetaData
   name::AbstractString
   fid::IOStream
   footer::XMLElement
   variable::Vector{String}
   "sorted cell IDs"
   cellid::Vector{UInt64}
   "ordered sequence index of raw cell IDs"
   cellIndex::Vector{Int64}
   time::Float64
   maxamr::Int64
   ncells::Vector{Int64}
   block_size::Vector{Int64}
   coordmin::Vector{Float64}
   coordmax::Vector{Float64}
   dcoord::Vector{Float64}
   populations::Vector{String}
   meshes::Dict{String, VMeshInfo}
end


function Base.show(io::IO, meta::MetaData)
   println(io, "filename = ", meta.name)
   println(io, "time = ", round(meta.time, digits=2))
   println(io, "dimension: $(ndims(meta))")
   println(io, "maximum AMR level: $(meta.maxamr)")
   println(io, "contains VDF: $(hasvdf(meta))")
   print(io, "variables: ")
   println(io, meta.variable)
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
   endian_offset = 8
   seek(fid, endian_offset)
   # Obtain the offset of the XML file
   offset = read(fid, UInt64)
   seek(fid, offset)
   footer = read(fid, String) |> parse_string |> root
end


"Return size and type information for the object."
function getObjInfo(fid, footer, name, tag, attr)

   arraysize, datasize, datatype, vectorsize, variable_offset = 0, 0, "", 0, 0
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

   !isFound && throw(ArgumentError("unknown variable $name"))

   seek(fid, variable_offset)

   if datatype == "float"
      T = datasize == 4 ? Float32 : Float64
   elseif datatype == "int"
      T = datasize == 4 ? Int32 : Int64
   elseif datatype == "uint"
      T = datasize == 4 ? UInt32 : UInt64
   end

   T, variable_offset, arraysize, datasize, vectorsize
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
    load(filename; verbose=false) -> MetaData

Return MetaData from a vlsv file.
"""
function load(filename::AbstractString; verbose=false)
   isfile(filename) || throw(ArgumentError("Cannot open \'$filename\': not a file"))
   fid = open(filename, "r")

   footer = getfooter(fid)

   meshName = "SpatialGrid"

   cellid = readvector(fid, footer, "CellID", "VARIABLE")

   cellIndex = sortperm(cellid)

   bbox = readmesh(fid, footer, meshName, "MESH_BBOX")

   nodeCoordsX = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_X")
   nodeCoordsY = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_Y")
   nodeCoordsZ = readmesh(fid, footer, meshName, "MESH_NODE_CRDS_Z")

   ncells = bbox[1:3]
   block_size = bbox[4:6]
   coordmin = [nodeCoordsX[begin], nodeCoordsY[begin], nodeCoordsZ[begin]]
   coordmax = [nodeCoordsX[end], nodeCoordsY[end], nodeCoordsZ[end]]

   dcoord = @. (coordmax - coordmin) / ncells

   meshes = Dict{String, VMeshInfo}()

   # Find all populations by the BLOCKIDS tag
   populations = String[]

   for varinfo in footer["BLOCKIDS"]

      if has_attribute(varinfo, "name")
         # VLSV 5.0 file with bounding box
         popname = attribute(varinfo, "name")

         bbox = readmesh(fid, footer, popname, "MESH_BBOX")

         nodeCoordsX = readmesh(fid, footer, popname, "MESH_NODE_CRDS_X")
         nodeCoordsY = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Y")
         nodeCoordsZ = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Z")
         vblocks = bbox[1:3]
         vblock_size = bbox[4:6]
         vmin = [nodeCoordsX[begin], nodeCoordsY[begin], nodeCoordsZ[begin]]
         vmax = [nodeCoordsX[end], nodeCoordsY[end], nodeCoordsZ[end]]
         dv = @. (vmax - vmin) / vblocks / vblock_size
      else
         popname = "avgs"

         if "vxblocks_ini" in attribute.(footer["PARAMETER"],"name")
            # In VLSV before 5.0 the mesh is defined with parameters.
            vblocks = @MVector zeros(Int, 3)
            vblocks[1] = readparameter(fid, footer, "vxblocks_ini")
            vblocks[2] = readparameter(fid, footer, "vyblocks_ini")
            vblocks[3] = readparameter(fid, footer, "vzblocks_ini")
            vblock_size = [4, 4, 4]
            vmin = [
               readparameter(fid, footer, "vxmin"),
               readparameter(fid, footer, "vymin"),
               readparameter(fid, footer, "vzmin") ]
            vmax = [
               readparameter(fid, footer, "vxmax"),
               readparameter(fid, footer, "vymax"),
               readparameter(fid, footer, "vzmax") ]
            dv = @. (vmax - vmin) / vblocks / vblock_size
         else
            # No velocity space info, e.g., file not written by Vlasiator
            vblocks = [0, 0, 0]
            vblock_size = [4, 4, 4]
            vmin, vmax, dv = zeros(3), zeros(3), ones(3)
         end
      end

      # Update list of active populations
      if popname ∉ populations
         push!(populations, popname)
      end

      # Create a new object for this population
      popVMesh = VMeshInfo(vblocks, vblock_size, vmin, vmax, dv)

      meshes[popname] = popVMesh

      verbose && @info "Found population $popname"
   end

   if hasname(footer, "PARAMETER", "t")
      timesim = readparameter(fid, footer, "t")
   elseif hasname(footer, "PARAMETER", "time") # Vlasiator 5.0+
      timesim = readparameter(fid, footer, "time")
   else
      timesim = Inf
   end

   # Obtain maximum refinement level
   ncell = prod(ncells)
   maxamr, cid = 0, ncell
   while cid < cellid[cellIndex[end]]
      maxamr += 1
      cid += ncell*8^maxamr
   end

   nVar = length(footer["VARIABLE"])
   vars = Vector{String}(undef, nVar)
   for i in 1:nVar
      @inbounds vars[i] = attribute(footer["VARIABLE"][i], "name")
   end

   # File IOstream is not closed for sake of data processing later.

   meta = MetaData(filename, fid, footer, vars, cellid[cellIndex], cellIndex, timesim,
      maxamr, ncells, block_size, coordmin, coordmax, dcoord, populations, meshes)
end


"""
    readvariablemeta(meta, var) -> VarInfo

Return VarInfo about `var` in the vlsv file linked to `meta`.
"""
function readvariablemeta(meta, var)

   unit, unitLaTeX, variableLaTeX, unitConversion = "", "", "", ""

   var = lowercase(var)

   # Get population and variable names from data array name
   if occursin("/", var)
      popname, varname = split(var, "/")
   else
      popname, varname = "pop", var
   end

   if hasvariable(meta, var) # For Vlasiator 5 vlsv files, metadata is included
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

   VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

"Return mesh related variable."
function readmesh(fid, footer, typeMesh, varMesh)
   T, _, arraysize, _, _ = getObjInfo(fid, footer, typeMesh, varMesh, "mesh")

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
   @unpack fid, footer, cellIndex = meta
   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
      return data
   end

   data = readvector(fid, footer, var, "VARIABLE")

   if startswith(var, "fg_") # fsgrid
      bbox = readmesh(fid, footer, "fsgrid", "MESH_BBOX")
      # Determine fsgrid domain decomposition
      nIORanks = readparameter(meta, "numWritingRanks")

      if ndims(data) > 1
         orderedData = zeros(Float32, size(data,1), bbox[1:3]...)
      else
         orderedData = zeros(Float32, bbox[1:3]...)
      end

      fgDecomposition = getDomainDecomposition(bbox[1:3], nIORanks)

      currentOffset = ones(Int, nIORanks+1)
      lsize = ones(Int, 3, nIORanks)
      lstart = similar(lsize)
      @inbounds @floop for i = 0:nIORanks-1
         x = i ÷ fgDecomposition[3] ÷ fgDecomposition[2]
         y = i ÷ fgDecomposition[3] % fgDecomposition[2]
         z = i % fgDecomposition[3]

         lsize[:,i+1] = calcLocalSize.(bbox[1:3], fgDecomposition, [x,y,z])
         lstart[:,i+1] = calcLocalStart.(bbox[1:3], fgDecomposition, [x,y,z])

         totalSize = prod(lsize[:,i+1])
         currentOffset[i+2] = currentOffset[i+1] + totalSize
      end

      @inbounds @floop for i = 1:nIORanks
         lend = lstart[:,i] + lsize[:,i] .- 1

         # Reorder data
         if ndims(data) > 1
            ldata = data[:,currentOffset[i]:currentOffset[i+1]-1]
            ldata = reshape(ldata, size(data,1), lsize[:,i]...)

            orderedData[:,lstart[1,i]:lend[1],lstart[2,i]:lend[2],lstart[3,i]:lend[3]] =
               ldata
         else
            ldata = data[currentOffset[i]:currentOffset[i+1]-1]
            ldata = reshape(ldata, lsize[:,i]...)

            orderedData[lstart[1,i]:lend[1],lstart[2,i]:lend[2],lstart[3,i]:lend[3]] = ldata
         end
      end
      data = dropdims(orderedData, dims=(findall(size(orderedData) .== 1)...,))
   elseif sorted # dccrg grid
      if ndims(data) == 1
         data = data[cellIndex]
      elseif ndims(data) == 2
         data = data[:, cellIndex]
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
   @unpack fid, footer = meta

   if isempty(ids)
      w = readvariable(meta, var, false)
      return [w]
   end

   T, offset, _, datasize, vectorsize = getObjInfo(fid, footer, var, "VARIABLE", "name")

   cellid = readvector(fid, footer, "CellID", "VARIABLE")
   rOffsets = [(findfirst(==(i), cellid)-1)*datasize*vectorsize for i in ids]

   v = Array{T}(undef, vectorsize, length(ids))

   for (i, r) in enumerate(rOffsets)
      seek(fid, offset + r)
      read!(fid, @view v[:,i])
   end

   return v
end

# Optimize decomposition of this grid over the given number of processors.
# Reference: fsgrid.hpp
function getDomainDecomposition(globalsize, nprocs)
   domainDecomp = @MVector [1, 1, 1]
   procBox = @MVector [0.0, 0.0, 0.0]
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

            v = 10*prod(procBox) + nyz + nzx + nxy

            if i * j * k == nprocs && v < minValue
               minValue = v
               domainDecomp[:] = [i, j, k]
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
    ndims(meta) -> Int

Return the dimension of VLSV data.
"""
Base.ndims(meta::MetaData) = count(>(1), meta.ncells)

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
   @unpack fid, footer = meta
   @unpack vblock_size = meta.meshes[pop]
   bsize = prod(vblock_size)

   cellsWithVDF = readvector(fid, footer, pop, "CELLSWITHBLOCKS")
   nblock_C = readvector(fid, footer, pop, "BLOCKSPERCELL")

   nblock_C_offsets = zeros(Int, length(cellsWithVDF))
   nblock_C_offsets[2:end] = cumsum(nblock_C[1:end-1])

   # Check if cells have vspace stored
   if cellid ∈ cellsWithVDF
      cellWithVDFIndex = findfirst(==(cellid), cellsWithVDF)
   else
      throw(ArgumentError("Cell ID $cellid does not store velocity distribution!"))
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
   variable_offset = parse(Int, content(varinfo))

   # Navigate to the correct position
   offset_f = offset * datasize + variable_offset

   seek(fid, offset_f)

   if datatype == "uint" && datasize == 4
      T = UInt32
   elseif datatype == "uint" && datasize == 8
      T = UInt64
   end

   blockIDs = Vector{T}(undef, nblocks)
   read!(fid, blockIDs)

   # Velocity cell IDs and corresponding avg distributions
   vcellids = zeros(Int, bsize*nblocks)
   vcellf = zeros(Tavg, bsize*nblocks)

   vcellid_local = [i + vblock_size[1]*j + vblock_size[1]*vblock_size[2]*k
      for i in 0:vblock_size[1]-1, j in 0:vblock_size[2]-1, k in 0:vblock_size[3]-1]

   @inbounds @floop for i in eachindex(blockIDs), j = 1:bsize
      vcellids[(i-1)*bsize+j] = vcellid_local[j] + bsize*blockIDs[i]
      vcellf[(i-1)*bsize+j] = data[j,i]
   end
   vcellids, vcellf
end


"""
    getvcellcoordinates(meta, vcellids, pop="proton")

Return velocity cells' coordinates of population `pop` and id `vcellids`.
"""
function getvcellcoordinates(meta, vcellids; pop="proton")
   @unpack vblocks, vblock_size, dv, vmin = meta.meshes[pop]

   bsize = prod(vblock_size)
   blockid = @. vcellids ÷ bsize
   # Get block coordinates
   blockIndX = @. blockid % vblocks[1]
   blockIndY = @. blockid ÷ vblocks[1] % vblocks[2]
   blockIndZ = @. blockid ÷ (vblocks[1] * vblocks[2])
   blockCoordX = @. blockIndX * dv[1] * vblock_size[1] + vmin[1]
   blockCoordY = @. blockIndY * dv[2] * vblock_size[2] + vmin[2]
   blockCoordZ = @. blockIndZ * dv[3] * vblock_size[3] + vmin[3]
   # Get cell indices
   cellids = @. vcellids % bsize
   cellidx = @. cellids % vblock_size[1]
   cellidy = @. cellids ÷ vblock_size[1] % vblock_size[2]
   cellidz = @. cellids ÷ (vblock_size[1] * vblock_size[2])
   # Get cell coordinates
   cellCoords = Matrix{Float32}(undef, 3, length(cellids))
   @inbounds @floop for i in eachindex(cellids)
      cellCoords[1,i] = blockCoordX[i] + (cellidx[i] + 0.5) * dv[1]
      cellCoords[2,i] = blockCoordY[i] + (cellidy[i] + 0.5) * dv[2]
      cellCoords[3,i] = blockCoordZ[i] + (cellidz[i] + 0.5) * dv[3]
   end
   cellCoords
end