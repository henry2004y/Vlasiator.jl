# VLSV reader in Julia
#
# Hongyang Zhou, hyzhou@umich.edu

include("vlsvvariables.jl")

using LightXML

export MetaData, VarInfo
export read_meta, read_variable, read_parameter, show_variables, has_variable,
       has_parameter, has_name, read_variable_select, read_variable_info,
       get_variable_derived, read_velocity_cells, get_velocity_cell_coordinates

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
Variable metadata from the vlsv footer, including the unit of the variable as a
regular string, the unit of the variable as a LaTeX-formatted string, the
description of the variable as a LaTeX-formatted string, and the 
conversion factor to SI units as a string.
"""
struct VarInfo
   unit::String
   unitLaTeX::LaTeXString
   variableLaTeX::LaTeXString
   unitConversion::String
end

"Meta data declaration."
struct MetaData
   name::AbstractString
   fid::IOStream
   footer::XMLElement
   cellid::Vector{UInt64}  # sorted cell IDs
   cellIndex::Vector{Int64}
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


function Base.show(io::IO, s::MetaData)
   println(io, "filename = ", s.name)
end

function Base.show(io::IO, s::VarInfo)
   println(io, "var in LaTeX: ", s.variableLaTeX)
   println(io, "unit: ", s.unit)
   println(io, "unit in LaTeX: ", s.unitLaTeX)
   println(io, "unit conversion = ", s.unitConversion)
end

"Return the xml footer of vlsv."
function read_xml_footer(fid)

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
function read_prep(fid, footer, name, tag, attr)

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

   return T::DataType, variable_offset, arraysize, datasize, vectorsize
end

"Return vector data from vlsv file."
function read_vector(fid, footer, name, tag)

   T, _, arraysize, _, vectorsize = read_prep(fid, footer, name, tag, "name")

   if vectorsize == 1
      w = Vector{T}(undef, arraysize)
   else
      w = Array{T,2}(undef, vectorsize, arraysize)
   end

   read!(fid, w)

   return w
end


"""
    read_meta(filename; verbose=false)

Return a struct of MetaData from vlsv file.
"""
function read_meta(filename::AbstractString; verbose=false)

   fid = open(filename, "r")

   footer = read_xml_footer(fid)

   meshName = "SpatialGrid"

   cellid = read_vector(fid, footer, "CellID", "VARIABLE")

   cellIndex = sortperm(cellid)

   bbox = read_mesh(fid, footer, meshName, "MESH_BBOX")

   nodeCoordsX = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_X")
   nodeCoordsY = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_Y")
   nodeCoordsZ = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_Z")
  
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

         bbox = read_mesh(fid, footer, popname, "MESH_BBOX")

         nodeCoordsX = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_X")   
         nodeCoordsY = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_Y")   
         nodeCoordsZ = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_Z")   
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
            vxblocks = read_parameter(fid, footer, "vxblocks_ini")
            vyblocks = read_parameter(fid, footer, "vyblocks_ini")
            vzblocks = read_parameter(fid, footer, "vzblocks_ini")
            vxblock_size = 4
            vyblock_size = 4
            vzblock_size = 4
            vxmin = read_parameter(fid, footer, "vxmin")
            vymin = read_parameter(fid, footer, "vymin")
            vzmin = read_parameter(fid, footer, "vzmin")
            vxmax = read_parameter(fid, footer, "vxmax")
            vymax = read_parameter(fid, footer, "vymax")
            vzmax = read_parameter(fid, footer, "vzmax")
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

   #close(fid) # Is it safe not to close it?

   meta = MetaData(filename, fid, footer, cellid[cellIndex], cellIndex,
      xcells, ycells, zcells, xblock_size, yblock_size, zblock_size,
      xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, meshes, populations)
end


"""
    read_variable_info(meta::MetaData, var)

Return a struct of VarInfo.           
"""
function read_variable_info(meta, var)

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

   if has_variable(meta.footer, var)
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

   return VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

"Return mesh related variable."
function read_mesh(fid, footer, typeMesh, varMesh)

   T, variable_offset, arraysize, datasize, vectorsize = 
      read_prep(fid, footer, typeMesh, varMesh, "mesh")

   w = Vector{T}(undef, arraysize)
   read!(fid, w)

   return w
end

"""
    read_variable(meta::MetaData, var, sorted=true)

Return variable value from the vlsv file. For DCCRG grid, the variables are
sorted by cell ID by default.
"""
function read_variable(meta, var, sorted=true)
   data = read_vector(meta.fid, meta.footer, var, "VARIABLE")
   
   if startswith(var, "fg_") # fsgrid
      bbox = read_mesh(meta.fid, meta.footer, "fsgrid", "MESH_BBOX")
      # Determine fsgrid domain decomposition
      nIORanks = read_parameter(meta, "numWritingRanks")

      if ndims(data) > 1
         orderedData = zeros(Float32, size(data,1), bbox[1:3]...)
      else
         orderedData = zeros(Float32, bbox[1:3]...)
      end

      currentOffset = 1
      fgDecomposition = getDomainDecomposition(bbox[1:3], nIORanks)

      for i = 0:nIORanks-1
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

            orderedData[:,
               lstart[1]:lend[1],lstart[2]:lend[2],lstart[3]:lend[3]] = ldata
         else
            ldata = data[currentOffset:currentOffset+totalSize-1]
            ldata = reshape(ldata, lsize...)

            orderedData[
               lstart[1]:lend[1],lstart[2]:lend[2],lstart[3]:lend[3]] = ldata
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

            value = 10*procBox[1]*procBox[2]*procBox[3] + nyz + nzx + nxy            

            if i * j * k == nprocs && value < minValue
               minValue = value
               domainDecomp[1:3] = [i, j, k]
            end
         end
      end
   end

   return domainDecomp
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

"Check if the VLSV file contains a variable."
has_variable(footer, var) = has_name(footer, "VARIABLE", var)

"""
    get_variable_derived(meta::MetaData, var)

Calculate a derived variable from basic quantities.
"""
function get_variable_derived(meta::MetaData, var)
   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
   else
      @error "Derived variable not known!"
   end
end

"""
    read_variable_select(meta::MetaData, var, idIn=UInt[])

Read a variable in a collection of cells.
"""
function read_variable_select(meta, var, idIn=UInt[])

   if isempty(idIn)
      w = read_variable(meta, var, false)
      return [w]
   end

   T, variable_offset, arraysize, datasize, vectorsize = 
      read_prep(meta.fid, meta.footer, var, "VARIABLE", "name")

   cellids = read_variable(meta, "CellID", false)
   rOffsets = [(findfirst(x->x==i, cellids)-1)*datasize*vectorsize for i in idIn]

   v = fill(T[], length(rOffsets))

   for (i, r) in enumerate(rOffsets)
      loc = variable_offset + r
      seek(meta.fid, loc)
   
      w = Vector{T}(undef, vectorsize)
      read!(meta.fid, w)
      v[i] = w
   end

   return v
end

"""
    read_parameter(meta::MetaData, param)

Return the parameter value from vlsv file.
"""
read_parameter(meta::MetaData, param) = read_parameter(meta.fid, meta.footer, param)

function read_parameter(fid, footer, param)
   
   T, _, _, _, _ = read_prep(fid, footer, param, "PARAMETER", "name")

   p = read(fid, T)
end

"""
    has_parameter(meta::MetaData, param)

Check if vlsv file contains a parameter.
"""
has_parameter(meta::MetaData, param) = has_name(meta.footer, "PARAMETER", param)

"Check if the XMLElement `elem` contains a `tag` with `name`."
function has_name(elem, tag, name)
   isFound = false
   
   for varinfo in elem[tag]
      attribute(varinfo, "name") == name && (isFound = true)
      isFound && break
   end
   
   return isFound
end

"Display all variables in the VLSV file."
function show_variables(meta::MetaData)
   nVar = length(meta.footer["VARIABLE"])
   vars = Vector{String}(undef, nVar)
   for i in 1:nVar
      vars[i] = attribute(meta.footer["VARIABLE"][i], "name")
   end
   vars
end

"""
    read_velocity_cells(meta, cellid; pop="proton")

Read velocity cells from a spatial cell of ID `cellid`, and return a map of
velocity cell ids and corresponding value.
"""
function read_velocity_cells(meta, cellid; pop="proton")

   vmesh = meta.meshes[pop]
   nblockx = vmesh.vxblock_size
   nblocky = vmesh.vyblock_size
   nblockz = vmesh.vzblock_size
   bsize = vmesh.vxblock_size * vmesh.vyblock_size * vmesh.vzblock_size

   cellsWithVDF = read_vector(meta.fid, meta.footer, pop, "CELLSWITHBLOCKS")
   nblock_C = read_vector(meta.fid, meta.footer, pop, "BLOCKSPERCELL")

   nblock_C_offsets = zeros(Int, length(cellsWithVDF))
   nblock_C_offsets[2:end] = cumsum(nblock_C[1:end-1])

   # Check that cells has vspace
   if cellid ∈ cellsWithVDF
      cellWithVDFIndex = findfirst(x->x==cellid, cellsWithVDF)
   else
      @error "The input cell does not have velocity distribution!"
   end

   # Navigate to the correct position
   offset = nblock_C_offsets[cellWithVDFIndex]
   nblocks = nblock_C[cellWithVDFIndex]

   # Read in avgs
   varinfo = meta.footer["BLOCKVARIABLE"][1]
   arraysize = parse(Int, attribute(varinfo, "arraysize"))
   datasize = parse(Int, attribute(varinfo, "datasize"))
   datatype = attribute(varinfo, "datatype")
   vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
   variable_offset = parse(Int, content(varinfo))

   # Navigate to the correct position
   offset_data = offset * vectorsize * datasize + variable_offset

   seek(meta.fid, offset_data)

   if datatype == "float" && datasize == 4
      Tavg = Float32
   elseif datatype == "float" && datasize == 8
      Tavg = Float64
   end

   data = Array{Tavg,2}(undef, vectorsize, nblocks)
   read!(meta.fid, data)

   # Read in block IDs
   varinfo = meta.footer["BLOCKIDS"][1]
   arraysize = parse(Int, attribute(varinfo, "arraysize"))
   datasize = parse(Int, attribute(varinfo, "datasize"))
   datatype = attribute(varinfo, "datatype")
   vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
   variable_offset = parse(Int, content(varinfo))

   # Navigate to the correct position
   offset_f = offset * vectorsize * datasize + variable_offset

   seek(meta.fid, offset_f)

   if datatype == "uint" && datasize == 4
      T = UInt32
   elseif datatype == "uint" && datasize == 8
      T = UInt64
   end

   blockIDs = Vector{T}(undef, nblocks*vectorsize)
   read!(meta.fid, blockIDs)

   # Velocity cell IDs and corresponding avg distributions
   vcellids = zeros(Int, bsize*nblocks)
   vcellf = zeros(Tavg, bsize*nblocks)

   vcellid_local = [i + nblockx*j + nblockx*nblocky*k
      for i in 0:nblockx-1, j in 0:nblocky-1, k in 0:nblockz-1]

   for i in 1:nblocks
      vblockid = blockIDs[i]
      for j = 1:bsize
         vcellids[(i-1)*bsize+j] = vcellid_local[j] + bsize*vblockid
         vcellf[(i-1)*bsize+j] = data[j,i]
      end
   end
   return vcellids, vcellf
end


"""
    get_velocity_cell_coordinates(meta, vcellids, pop="proton")

Returns velocity cells' coordinates of population `pop` and id `vcellids`.
"""
function get_velocity_cell_coordinates(meta, vcellids; pop="proton")

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
   #cellCoords = [blockCoordX + (cellidx + 0.5) * vmesh.dvx,
   #              blockCoordY + (cellidy + 0.5) * vmesh.dvy,
   #              blockCoordZ + (cellidz + 0.5) * vmesh.dvz]

   cellCoords = Matrix{Float32}(undef, 3, length(cellids))
   for i = 1:length(cellids)
      cellCoords[1,i] = blockCoordX[i] + (cellidx[i] + 0.5) * vmesh.dvx
      cellCoords[2,i] = blockCoordY[i] + (cellidy[i] + 0.5) * vmesh.dvy
      cellCoords[3,i] = blockCoordZ[i] + (cellidz[i] + 0.5) * vmesh.dvz
   end
   cellCoords
end