# VLSV reader in Julia

include("vlsvvariables.jl")

"Velocity mesh information."
struct VMeshInfo
   "number of velocity blocks"
   vblocks::NTuple{3, Int64}
   vblock_size::NTuple{3, Int64}
   vmin::NTuple{3, Float64}
   vmax::NTuple{3, Float64}
   dv::NTuple{3, Float64}
end

"Variable information from the VLSV footer."
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
struct MetaVLSV
   dir::String
   name::String
   fid::IOStream
   footer::EzXML.Node
   variable::Vector{String}
   "mapping of unsorted cell ID to ordering"
   celldict::Dict{UInt64, Int64}
   "ordered sequence indexes of raw cell IDs"
   cellindex::Vector{Int64}
   time::Float64
   maxamr::Int64
   hasvdf::Bool
   ncells::NTuple{3, Int64}
   block_size::NTuple{3, Int64}
   coordmin::NTuple{3, Float64}
   coordmax::NTuple{3, Float64}
   dcoord::NTuple{3, Float64}
   species::Vector{String}
   meshes::Dict{String, VMeshInfo}
end


function Base.show(io::IO, meta::MetaVLSV)
   println(io, "File: ", meta.name)
   println(io, "Time: ", round(meta.time, digits=2))
   println(io, "Dimension: ", ndims(meta))
   println(io, "Maximum AMR level: ", meta.maxamr)
   println(io, "Contains VDF: ", meta.hasvdf)
   println(io, "Variables: ", meta.variable)
end

function Base.show(io::IO, s::VarInfo)
   println(io, "Variable in LaTeX: ", s.variableLaTeX)
   println(io, "Unit: ", s.unit)
   println(io, "Unit in LaTeX: ", s.unitLaTeX)
   println(io, "Unit conversion: ", s.unitConversion)
end

function Base.show(io::IO, vmesh::VMeshInfo)
   println(io, "vblocks: ", vmesh.vblocks)
   println(io, "vblock size: ", vmesh.vblock_size)
   foreach((vmin,dv,vmax,comp) -> println(io, "v$comp range: ", vmin, ":", dv, ":", vmax),
      vmesh.vmin, vmesh.dv, vmesh.vmax, 'x':'z')
end

"Return the XML footer of opened VLSV file."
@inline function getfooter(fid::IOStream)
   # First 8 bytes indicate big-endian or else
   endian_offset = 8
   seek(fid, endian_offset)
   # Obtain the offset of the XML footer
   offset = read(fid, UInt64)
   seek(fid, offset)
   footer = read(fid, String) |> parsexml |> root
end


"Return size and type information for the inquired parameter/data with `name`."
function getObjInfo(footer::EzXML.Node, name::String, tag::String, attr::String)
   local arraysize, datasize, datatype, vectorsize, variable_offset
   isFound = false

   for var in findall("//$tag", footer)
      if var[attr] == String(name)
         arraysize = parse(Int, var["arraysize"])
         datasize = parse(Int, var["datasize"])
         datatype = var["datatype"]
         vectorsize = parse(Int, var["vectorsize"])
         variable_offset = parse(Int, nodecontent(var))
         isFound = true
         break
      end
   end

   isFound || throw(ArgumentError("unknown variable $name"))

   T::Type =
      if datatype == "float"
         datasize == 4 ? Float32 : Float64
      elseif datatype == "int"
         datasize == 4 ? Int32 : Int64
      elseif datatype == "uint"
         datasize == 4 ? UInt32 : UInt64
      end

   T, variable_offset, arraysize, datasize, vectorsize
end

"Return vector of `name` from the VLSV file with `footer` associated with stream `fid`."
function readvector(fid::IOStream, footer::EzXML.Node, name::String, tag::String)
   T, offset, asize, dsize, vsize = getObjInfo(footer, name, tag, "name")

   if Sys.total_memory() > 8*asize*vsize*dsize
      w = vsize == 1 ?
         Vector{T}(undef, asize) :
         Array{T,2}(undef, vsize, asize)
      seek(fid, offset)
      read!(fid, w)
   else
      @warn "Large array detected. Using memory-mapped I/O!" maxlog=1
      a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
      w = vsize == 1 ?
         reinterpret(T, a) :
         reshape(reinterpret(T, a), vsize, asize)
   end

   w
end

"""
    load(file::AbstractString)) -> MetaVLSV

Generate a MetaVLSV object from `file` of VLSV format.
"""
function load(file::AbstractString)
   fid = open(file, "r")

   footer = getfooter(fid)

   cellid = getcellid(fid, footer)

   cellindex = sortperm(cellid)

   bbox = readmesh(fid, footer, "SpatialGrid", "MESH_BBOX")::Vector{UInt}

   nodeX = readmesh(fid, footer, "SpatialGrid", "MESH_NODE_CRDS_X")::Vector{Float64}
   nodeY = readmesh(fid, footer, "SpatialGrid", "MESH_NODE_CRDS_Y")::Vector{Float64}
   nodeZ = readmesh(fid, footer, "SpatialGrid", "MESH_NODE_CRDS_Z")::Vector{Float64}

   @inbounds ncells = (bbox[1], bbox[2], bbox[3])
   @inbounds block_size = (bbox[4], bbox[5], bbox[6])
   @inbounds coordmin = (nodeX[begin], nodeY[begin], nodeZ[begin])
   @inbounds coordmax = (nodeX[end], nodeY[end], nodeZ[end])

   dcoord = ntuple(i -> (coordmax[i] - coordmin[i]) / ncells[i], Val(3))

   meshes = Dict{String, VMeshInfo}()

   # Find all species by the BLOCKIDS tag
   species = String[]

   vblocks = [0, 0, 0]
   vblock_size = [4, 4, 4]
   vmin = [0.0, 0.0, 0.0]
   vmax = [1.0, 1.0, 1.0]

   for varinfo in findall("//BLOCKIDS", footer)

      if haskey(varinfo, "name")
         # VLSV 5.0 file with bounding box
         popname = varinfo["name"]

         vbox = readmesh(fid, footer, popname, "MESH_BBOX")::Vector{UInt}
         nodeX = readmesh(fid, footer, popname, "MESH_NODE_CRDS_X")::Vector{Float64}
         nodeY = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Y")::Vector{Float64}
         nodeZ = readmesh(fid, footer, popname, "MESH_NODE_CRDS_Z")::Vector{Float64}

         vblocks[1], vblocks[2], vblocks[3] = vbox[1], vbox[2], vbox[3]
         vblock_size[1], vblock_size[2], vblock_size[3] = vbox[4], vbox[5], vbox[6]
         vmin[1], vmin[2], vmin[3] = nodeX[begin], nodeY[begin], nodeZ[begin]
         vmax[1], vmax[2], vmax[3] = nodeX[end], nodeY[end], nodeZ[end]
         dv = ntuple(i -> (vmax[i] - vmin[i]) / vblocks[i] / vblock_size[i], Val(3))
      else
         popname = "avgs"
         # In VLSV before 5.0 the mesh is defined with parameters.
         if "vxblocks_ini" in getindex.(findall("//PARAMETER", footer), "name")
            vblocks_str = ("vxblocks_ini", "vyblocks_ini", "vzblocks_ini")
            vmin_str = ("vxmin", "vymin", "vzmin")
            vmax_str = ("vxmax", "vymax", "vzmax")
            vblocks .= [readparameter(fid, footer, vblocks_str[i]) for i in 1:3]
            vmin .= [readparameter(fid, footer, vmin_str[i]) for i in 1:3]
            vmax .= [readparameter(fid, footer, vmax_str[i]) for i in 1:3]
            dv = ntuple(i -> (vmax[i] - vmin[i]) / vblocks[i] / vblock_size[i], Val(3))
         else
            # No velocity space info, e.g., file not written by Vlasiator
            dv = (1.0, 1.0, 1.0)
         end
      end

      # Update list of active species
      if popname ∉ species
         push!(species, popname)
      end

      # Create a new object for this population
      popVMesh = VMeshInfo(
         (vblocks[1], vblocks[2], vblocks[3]),
         (vblock_size[1], vblock_size[2], vblock_size[3]),
         (vmin[1], vmin[2], vmin[3]),
         (vmax[1], vmax[2], vmax[3]),
         dv)

      meshes[popname] = popVMesh
   end

   if hasname(footer, "PARAMETER", "time") # Vlasiator 5.0+
      timesim = readparameter(fid, footer, "time")::Float64
   elseif hasname(footer, "PARAMETER", "t")
      timesim = readparameter(fid, footer, "t")::Float64
   else
      timesim = Inf
   end

   celldict = Dict(cellid[i] => i for i in eachindex(cellid))

   # Obtain maximum refinement level
   maxamr = getmaxrefinement(cellid, ncells)

   varinfo = findall("//VARIABLE", footer)
   vars = [info["name"] for info in varinfo]

   hasvdf = let
      vcells = readmesh(fid, footer, "SpatialGrid", "CELLSWITHBLOCKS")::Vector{UInt}
      !isempty(vcells)
   end

   # File IOstream is not closed for sake of data processing later.

   meta = MetaVLSV(splitdir(file)..., fid, footer, vars, celldict, cellindex,
      timesim, maxamr, hasvdf, ncells, block_size, coordmin, coordmax,
      dcoord, species, meshes)
end

# Allow `do ... end` syntax.
function load(f::Function, file::AbstractString)
   meta = load(file)
   try
      f(meta)
   finally
      close(meta.fid)
   end
end

function getmaxrefinement(cellid::AbstractArray{UInt64, 1}, ncells::NTuple{3, UInt64})

   ncell = prod(ncells)
   maxamr, cid = 0, ncell
   while @inbounds cid < maximum(cellid)
      maxamr += 1
      cid += ncell*8^maxamr
   end

   maxamr
end

"""
    readvariablemeta(meta, var) -> VarInfo

Return VarInfo about `var` in the VLSV file associated with `meta`.
"""
function readvariablemeta(meta::MetaVLSV, var::String)
   varSym = isa(var, AbstractString) ? Symbol(var) : var

   unit, unitLaTeX, variableLaTeX, unitConversion = "", "", "", ""

   if varSym in keys(units_predefined)
      unit, variableLaTeX, unitLaTeX = units_predefined[varSym]
   elseif hasvariable(meta, var) # For Vlasiator 5 files, MetaVLSV is included
      for varinfo in findall("//VARIABLE", meta.footer)
         if varinfo["name"] == var
            haskey(varinfo, "unit") || break
            unit = varinfo["unit"]
            unitLaTeX = varinfo["unitLaTeX"]
            variableLaTeX = varinfo["variableLaTeX"]
            unitConversion = varinfo["unitConversion"]
         end
      end
   end

   VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

"Return mesh related variable."
function readmesh(fid::IOStream, footer::EzXML.Node, typeMesh::String, varMesh::String)
   T, offset, arraysize, _, _ = getObjInfo(footer, typeMesh, varMesh, "mesh")

   w = Vector{T}(undef, arraysize)
   seek(fid, offset)
   read!(fid, w)

   w
end

"""
    readvariable(meta::MetaVLSV, var::String, sorted::Bool=true) -> Array

Return variable value of `var` from the VLSV file associated with `meta`. By default
`sorted=true`, which means that for DCCRG grid the variables are sorted by cell ID.
"""
function readvariable(meta::MetaVLSV, var::String, sorted::Bool=true)
   (;fid, footer, cellindex) = meta
   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta)::Array
      return v
   end

   raw = readvector(fid, footer, var, "VARIABLE")

   if startswith(var, "fg_") # fsgrid
      bbox = readmesh(fid, footer, "fsgrid", "MESH_BBOX")::Vector{Int}
      # Determine fsgrid domain decomposition
      nIORanks = readparameter(meta, "numWritingRanks")::Int32

      dataOrdered =
         if ndims(raw) > 1
            @inbounds Array{Float32}(undef, size(raw,1), bbox[1], bbox[2], bbox[3])
         else
            @inbounds Array{Float32}(undef, bbox[1], bbox[2], bbox[3])
         end

      @inbounds fgDecomposition = @views getDomainDecomposition(bbox[1:3], nIORanks)

      _fillFGordered!(dataOrdered, raw, fgDecomposition, nIORanks, bbox)

      v = dropdims(dataOrdered, dims=(findall(size(dataOrdered) .== 1)...,))
   elseif sorted # dccrg grid
      @inbounds v = ndims(raw) == 1 ? raw[cellindex] : raw[:,cellindex]
      if eltype(v) == Float64; v = Float32.(v); end
   else
      v = raw
   end

   return v::Array
end

"""
    readvariable(meta::MetaVLSV, var::String, ids::Vector{<:Integer}) -> Array

Read variable `var` in a collection of cells `ids` associated with `meta`. if `ids` is
empty, return the whole sorted array of `var`.
"""
function readvariable(meta::MetaVLSV, var::String,
   ids::Union{Vector{<:Integer}, UnitRange{Int64}})

   startswith(var, "fg_") && error("Currently does not support reading fsgrid!")
   (;fid, footer, celldict) = meta

   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta, ids)
      return v::Array
   end

   if isempty(ids)
      v = readvariable(meta, var)
      return v::Array
   else
      T, offset, asize, dsize, vsize = getObjInfo(footer, var, "VARIABLE", "name")

      v = Array{T}(undef, vsize, length(ids))

      w = let
         a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
         reshape(reinterpret(T, a), vsize, asize)
      end

      id_ = [celldict[id] for id in ids]

      _fillv!(v, w, id_, vsize)

      if T == Float64
         v = Float32.(v)
      end

      return v::Array
   end
end

"""
    readvariable(meta::MetaVLSV, var::String, cid::Integer) -> Array

Read variable `var` in cell `cid` associated with `meta`.
"""
function readvariable(meta::MetaVLSV, var::String, cid::Integer)
   startswith(var, "fg_") && error("Currently does not support reading fsgrid!")
   (;fid, footer, celldict) = meta

   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta, cid)::Array
      return v
   end

   T, offset, _, dsize, vsize = getObjInfo(footer, var, "VARIABLE", "name")

   v = vsize == 1 ? Vector{T}(undef, 1) : Array{T,2}(undef, vsize, 1)
   seek(fid, offset + (celldict[cid]-1)*vsize*dsize)
   read!(fid, v)

   if T == Float64
      v = Float32.(v)
   end

   return v::Array
end


@inline function _fillv!(v, w, id_, vsize)
   for i in eachindex(id_), iv = 1:vsize
      @inbounds v[iv,i] = w[iv,id_[i]]
   end
end

function _fillFGordered!(dataOrdered, raw, fgDecomposition, nIORanks, bbox)
   offsetnow = 1

   @inbounds @views for i = 1:nIORanks
      xyz = (
         (i - 1) ÷ fgDecomposition[3] ÷ fgDecomposition[2],
         (i - 1) ÷ fgDecomposition[3] % fgDecomposition[2],
         (i - 1) % fgDecomposition[3] )

      lsize = ntuple(i -> calcLocalSize(bbox[i], fgDecomposition[i], xyz[i]), Val(3))
      lstart = ntuple(i -> calcLocalStart(bbox[i], fgDecomposition[i], xyz[i]), Val(3))
      offsetnext = offsetnow + prod(lsize)
      lend = @. lstart + lsize - 1
      lrange = map((x,y)->x:y, lstart, lend)
      # Reorder data
      if ndims(raw) > 1
         ldata = reshape(raw[:,offsetnow:offsetnext-1], size(raw,1), lsize...)

         dataOrdered[:,lrange...] = ldata
      else
         ldata = reshape(raw[offsetnow:offsetnext-1], lsize...)

         dataOrdered[lrange...] = ldata
      end
      offsetnow = offsetnext
   end

   return
end

@inline @Base.propagate_inbounds Base.getindex(meta::MetaVLSV, key::String) =
   readvariable(meta, key)

@inline function getcellid(fid::IOStream, footer::EzXML.Node)
   _, offset, asize, _, _ = getObjInfo(footer, "CellID", "VARIABLE", "name")
   a = mmap(fid, Vector{UInt8}, 8*asize, offset)
   cellid = reinterpret(UInt64, a)
end

"""
    extractsat(files::Vector{String}, var::String, cid::Integer)

Extract `var` at a fixed cell ID `cid` from `files`. This assumes that `files` come from the
same grid structure.
"""
function extractsat(files::AbstractVector{String}, var::String, cid::Integer)
   v = open(files[1], "r") do fid
      footer = getfooter(fid)
      T, _, _, _, vsize = getObjInfo(footer, var, "VARIABLE", "name")
      Array{T,2}(undef, vsize, length(files))
   end

   for i in eachindex(files)
      fid = open(files[i], "r")
      footer = getfooter(fid)
      cellid = getcellid(fid, footer)
      c_ = findfirst(isequal(cid), cellid)
      _, offset, _, dsize, vsize = getObjInfo(footer, var, "VARIABLE", "name")
      seek(fid, offset + (c_ - 1)*vsize*dsize)
      @views read!(fid, v[:,i])
   end

   v
end

"""
    getdata2d(meta::MetaVLSV, var::String)

Return 2d scalar/vector data. Nonpublic since it won't work with DCCRG AMR.
"""
function getdata2d(meta::MetaVLSV, var::String)
   ndims(meta) == 2 || @error "2D outputs required."
   sizes = filter(!=(1), meta.ncells)
   data = readvariable(meta, var)
   data = ndims(data) == 1 || size(data, 1) == 1 ?
      reshape(data, sizes[1], sizes[2]) :
      reshape(data, 3, sizes[1], sizes[2]) # assumes 3-vector, may not work in general
end

"File size in bytes."
@inline Base.size(meta::MetaVLSV) = filesize(joinpath(meta.dir, meta.name))

"""
    getDomainDecomposition(globalsize, nprocs)

Obtain decomposition of this grid over the given number of processors.
Reference: fsgrid.hpp
"""
function getDomainDecomposition(globalsize, nprocs)
   domainDecomp = (1, 1, 1)
   minValue = typemax(Float64)

   @inbounds for i = 1:min(nprocs, globalsize[1])
      iBox = max(globalsize[1]/i, 1.0)

      for j = 1:min(nprocs, globalsize[2])
         i * j > nprocs && break

         jBox = max(globalsize[2]/j, 1.0)

         for k = 1:min(nprocs, globalsize[2])
            i * j * k > nprocs && continue

            kBox = max(globalsize[3]/k, 1.0)

            nyz = ifelse(i > 1, jBox * kBox, 0.0)
            nzx = ifelse(j > 1, kBox * iBox, 0.0)
            nxy = ifelse(k > 1, iBox * jBox, 0.0)

            v = 10*iBox*jBox*kBox + nyz + nzx + nxy

            if i * j * k == nprocs && v < minValue
               minValue = v
               domainDecomp = (i, j, k)
            end
         end
      end
   end

   domainDecomp
end

function calcLocalStart(globalCells, nprocs, lcells)
   ncells = globalCells ÷ nprocs
   remainder = globalCells % nprocs
   lstart = ifelse(lcells < remainder, lcells*(ncells+1) + 1, lcells*ncells + remainder + 1)
end

function calcLocalSize(globalCells, nprocs, lcells)
   ncells = globalCells ÷ nprocs
   remainder = globalCells % nprocs
   lsize = ifelse(lcells < remainder, ncells + 1, ncells)
end

"""
    hasvariable(meta::MetaVLSV, var::String) -> Bool

Check if the VLSV file associated with `meta` contains a variable `var`.
"""
hasvariable(meta::MetaVLSV, var::String) = hasname(meta.footer, "VARIABLE", var)

"""
    readparameter(meta::MetaVLSV, param::String)

Return the parameter value from the VLSV file associated with `meta`.
"""
readparameter(meta::MetaVLSV, param::String) = readparameter(meta.fid, meta.footer, param)

function readparameter(fid::IOStream, footer::EzXML.Node, param::String)
   T, offset, _, _, _ = getObjInfo(footer, param, "PARAMETER", "name")
   seek(fid, offset)
   p = read(fid, T)
end

"""
    hasparameter(meta::MetaVLSV, param::String) -> Bool

Check if the VLSV file contains a certain parameter `param`.
"""
hasparameter(meta::MetaVLSV, param::String) = hasname(meta.footer, "PARAMETER", param)

"Check if the XML `element` contains a `tag` with `name`."
function hasname(element::EzXML.Node, tag::String, name::String)
   isFound = false

   for var in findall("//$tag", element)
      var["name"] == name && (isFound = true)
      isFound && break
   end

   isFound
end

"""
    ndims(meta::MetaVLSV) -> Int

Return the simulation spatial dimension of VLSV data.
"""
Base.ndims(meta::MetaVLSV) = count(>(1), meta.ncells)

"Determine whether `meta` is not yet closed."
Base.isopen(meta::MetaVLSV) = isopen(meta.fid)

"""
    readvcells(meta::MetaVLSV, cid::Integer; species="proton") -> vcellids, vcellf

Read velocity cells of `species` from a spatial cell of ID `cid` associated with `meta`, and
return a map of velocity cell ids `vcellids` and corresponding value `vcellf`.
"""
function readvcells(meta::MetaVLSV, cid::Integer; species::String="proton")
   (;fid, footer) = meta
   (;vblock_size) = meta.meshes[species]
   bsize = prod(vblock_size)

   local offset::Int, nblocks::Int
   let cellsWithVDF = readvector(fid, footer, species, "CELLSWITHBLOCKS"),
      nblock_C = readvector(fid, footer, species, "BLOCKSPERCELL")::Vector{UInt32}

      cellWithVDFIndex = findfirst(==(cid), cellsWithVDF)
      if !isnothing(cellWithVDFIndex)
         @inbounds nblocks = Int(nblock_C[cellWithVDFIndex])
         if nblocks == 0
            throw(ArgumentError("Cell ID $cid does not store velocity distribution!"))
         end
      else
         throw(ArgumentError("Cell ID $cid does not store velocity distribution!"))
      end
      # Offset position to vcell storage
      offset = sum(@view nblock_C[1:cellWithVDFIndex-1])
   end

   local datasize, datatype, vectorsize, variable_offset
   # Read in avgs
   let varinfo = findfirst("//BLOCKVARIABLE", footer)
      datasize = parse(Int, varinfo["datasize"])
      datatype = varinfo["datatype"]
      datatype == "float" || @error "VDFs must be floating numbers!"
      vectorsize = parse(Int, varinfo["vectorsize"])
      variable_offset = parse(Int, nodecontent(varinfo))
   end

   data = let
      Tavg = datasize == 4 ? Float32 : Float64
      a = mmap(fid, Vector{UInt8}, datasize*vectorsize*nblocks,
         offset*vectorsize*datasize + variable_offset)
      reshape(reinterpret(Tavg, a), vectorsize, nblocks)
   end

   # Read in block IDs
   let varinfo = findfirst("//BLOCKIDS", footer)
      datasize = parse(Int, varinfo["datasize"])
      datatype = varinfo["datatype"]
      datatype == "uint" || @error "Block ID must be unsigned integer!"
      variable_offset = parse(Int, nodecontent(varinfo))
   end

   blockIDs = let T = datasize == 4 ? UInt32 : UInt64
      Vector{T}(undef, nblocks)
   end
   seek(fid, offset * datasize + variable_offset)
   read!(fid, blockIDs)

   # Velocity cell IDs and distributions (ordered by blocks)
   vcellids = Vector{UInt32}(undef, bsize*nblocks)
   vcellf = Vector{Float32}(undef, bsize*nblocks)

   vcellid_local = @inbounds [i + vblock_size[1]*j + vblock_size[1]*vblock_size[2]*k
      for i in 0:vblock_size[1]-1, j in 0:vblock_size[2]-1, k in 0:vblock_size[3]-1]

   _fillvcell!(vcellids, vcellf, vcellid_local, data, blockIDs, bsize)

   vcellids, vcellf
end

@inline function _fillvcell!(vcellids, vcellf, vcellid_local, data, blockIDs, bsize)
   @inbounds for i in eachindex(blockIDs), j = 1:bsize
      index_ = (i-1)*bsize+j
      vcellids[index_] = vcellid_local[j] + bsize*blockIDs[i]
      vcellf[index_] = data[j,i]
   end
end