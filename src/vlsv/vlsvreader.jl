# VLSV reader in Julia

include("vlsvvariables.jl")

const NodeVector = SubArray{EzXML.Node, 1, Vector{EzXML.Node}, Tuple{UnitRange{Int64}},
   true}

"Velocity mesh information."
struct VMeshInfo
   "number of velocity blocks"
   vblocks::NTuple{3, Int}
   vblock_size::NTuple{3, Int}
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

struct NodeVLSV
   var::NodeVector
   param::NodeVector
   cellwithVDF::NodeVector
   cellblocks::NodeVector
   blockvar::NodeVector
   blockid::NodeVector
end

"VLSV meta data."
struct MetaVLSV
   dir::String
   name::String
   fid::IOStream
   nodeVLSV::NodeVLSV
   variable::Vector{String}
   "mapping of unsorted cell ID to ordering"
   celldict::Dict{UInt, Int}
   "ordered sequence indexes of raw cell IDs"
   cellindex::Vector{Int}
   time::Float64
   maxamr::Int
   hasvdf::Bool
   ncells::NTuple{3, Int}
   block_size::NTuple{3, Int}
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
   # Obtain the offset of the XML
   offset = read(fid, UInt)
   seek(fid, offset)
   footer = read(fid, String) |> parsexml |> root
end


@inline function getdatatype(datatype::String, datasize::Int)
   T::Type =
      if datatype == "float"
         datasize == 4 ? Float32 : Float64
      elseif datatype == "int"
         datasize == 4 ? Int32 : Int
      elseif datatype == "uint"
         datasize == 4 ? UInt32 : UInt
      else
         throw(ArgumentError("unknown type $datatype"))
      end
end

function getvarinfo(nodevar::AbstractVector{EzXML.Node}, name::String)
   local arraysize, datasize, datatype, vectorsize, offset
   isFound = false

   for var in nodevar
      if var["name"] == name
         arraysize = parse(Int, var["arraysize"])
         datasize = parse(Int, var["datasize"])
         datatype = var["datatype"]
         vectorsize = parse(Int, var["vectorsize"])
         offset = parse(Int, nodecontent(var))
         isFound = true
         break
      end
   end

   isFound || throw(ArgumentError("unknown variable $name"))

   T = getdatatype(datatype, datasize)

   T, offset, arraysize, datasize, vectorsize
end

@inline function getparaminfo(nodeparam::NodeVector, name::String)
   local datasize, datatype, offset
   isFound = false

   for p in nodeparam
      if p["name"] == name
         datasize = parse(Int, p["datasize"])
         datatype = p["datatype"]
         offset = parse(Int, nodecontent(p))
         isFound = true
         break
      end
   end

   isFound || throw(ArgumentError("unknown variable $name"))

   T = getdatatype(datatype, datasize)

   T, offset
end

"General inquiry of element `tag` with `name` and `attr`."
function getObjInfo(footer::EzXML.Node, name::String, tag::String, attr::String)
   local arraysize, datasize, datatype, vectorsize, offset
   isFound = false

   for var in findall("//$tag", footer)
      if var[attr] == name
         arraysize = parse(Int, var["arraysize"])
         datasize = parse(Int, var["datasize"])
         datatype = var["datatype"]
         vectorsize = parse(Int, var["vectorsize"])
         offset = parse(Int, nodecontent(var))
         isFound = true
         break
      end
   end

   isFound || throw(ArgumentError("unknown variable $name"))

   T = getdatatype(datatype, datasize)

   T, offset, arraysize, datasize, vectorsize
end

"Return vector of `name` from the VLSV file associated with stream `fid`."
@inline function readvector(fid::IOStream, ::Type{T}, offset::Int, asize::Int, dsize::Int,
   usemmap::Bool=false) where T
   if !usemmap
      w = Vector{T}(undef, asize)
      seek(fid, offset)
      read!(fid, w)
   else
      a = mmap(fid, Vector{UInt8}, dsize*asize, offset)
      w = reinterpret(T, a)
   end

   w::AbstractVector{T}
end

@inline function readarray(fid::IOStream, ::Type{T}, offset::Int, asize::Int, dsize::Int,
   vsize::Int, usemmap::Bool=false) where T
   if !usemmap
      w = Array{T,2}(undef, vsize, asize)
      seek(fid, offset)
      read!(fid, w)
   else
      a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
      w = reshape(reinterpret(T, a), vsize, asize)
   end

   w::AbstractArray{T, 2}
end

"""
    load(file::AbstractString)) -> MetaVLSV

Generate a MetaVLSV object from `file` of VLSV format.
"""
function load(file::AbstractString)
   fid = open(file, "r")

   footer = getfooter(fid)

   nodes = elements(footer)

   ibegin_, iend_ = zeros(Int, 6), zeros(Int, 6)

   for i in eachindex(nodes)
      if nodes[i].name == "VARIABLE"
         if ibegin_[1] == 0 ibegin_[1] = i end
         iend_[1] = i
      elseif nodes[i].name == "PARAMETER"
         if ibegin_[2] == 0 ibegin_[2] = i end
         iend_[2] = i
      elseif nodes[i].name == "CELLSWITHBLOCKS"
         if ibegin_[3] == 0 ibegin_[3] = i end
         iend_[3] = i
      elseif nodes[i].name == "BLOCKSPERCELL"
         if ibegin_[4] == 0 ibegin_[4] = i end
         iend_[4] = i
      elseif nodes[i].name == "BLOCKVARIABLE"
         if ibegin_[5] == 0 ibegin_[5] = i end
         iend_[5] = i
      elseif nodes[i].name == "BLOCKIDS"
         if ibegin_[6] == 0 ibegin_[6] = i end
         iend_[6] = i
      end
   end

   @views begin
      nodevar = nodes[ibegin_[1]:iend_[1]]
      nodeparam = nodes[ibegin_[2]:iend_[2]]
      nodecellwithVDF = nodes[ibegin_[3]:iend_[3]]
      nodecellblocks = nodes[ibegin_[4]:iend_[4]]
      nodeblockvar = nodes[ibegin_[5]:iend_[5]]
      nodeblockid = nodes[ibegin_[6]:iend_[6]]
   end

   n = NodeVLSV(nodevar, nodeparam, nodecellwithVDF, nodecellblocks, nodeblockvar,
      nodeblockid)

   cellid = getcellid(fid, n.var)

   cellindex = sortperm(cellid)

   ncells, block_size, coordmin, coordmax = readmesh(fid, footer)

   dcoord = ntuple(i -> (coordmax[i] - coordmin[i]) / ncells[i], Val(3))

   meshes = Dict{String, VMeshInfo}()

   # Find all species by the BLOCKIDS tag
   species = String[]

   vblocks = [0, 0, 0]
   vblock_size = [4, 4, 4]
   vmin = [0.0, 0.0, 0.0]
   vmax = [1.0, 1.0, 1.0]

   for node in n.blockid
      if haskey(node, "name")
         # VLSV 5.0 file with bounding box
         popname = node["name"]

         vbox, nodeX, nodeY, nodeZ = readvmesh(fid, footer, popname)

         vblocks[1], vblocks[2], vblocks[3] = vbox[1], vbox[2], vbox[3]
         vblock_size[1], vblock_size[2], vblock_size[3] = vbox[4], vbox[5], vbox[6]
         vmin[1], vmin[2], vmin[3] = nodeX[begin], nodeY[begin], nodeZ[begin]
         vmax[1], vmax[2], vmax[3] = nodeX[end], nodeY[end], nodeZ[end]
         dv = ntuple(i -> (vmax[i] - vmin[i]) / vblocks[i] / vblock_size[i], Val(3))
      else
         popname = "avgs"
         # In VLSV before 5.0 the mesh is defined with parameters.
         if "vxblocks_ini" in getindex.(n.param, "name")
            vblocks_str = ("vxblocks_ini", "vyblocks_ini", "vzblocks_ini")
            vmin_str = ("vxmin", "vymin", "vzmin")
            vmax_str = ("vxmax", "vymax", "vzmax")
            vblocks .= [readparameter(fid, n.param, vblocks_str[i]) for i in 1:3]
            vmin .= [readparameter(fid, n.param, vmin_str[i]) for i in 1:3]
            vmax .= [readparameter(fid, n.param, vmax_str[i]) for i in 1:3]
            dv = ntuple(i -> (vmax[i] - vmin[i]) / vblocks[i] / vblock_size[i], Val(3))
         else
            error("File not written by Vlasiator!")
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

   if hasname(n.param, "time") # Vlasiator 5.0+
      timesim = readparameter(fid, n.param, "time")::Float64
   elseif hasname(n.param, "t")
      timesim = readparameter(fid, n.param, "t")::Float64
   else
      timesim = Inf
   end

   celldict = Dict(cellid[i] => i for i in eachindex(cellid))

   # Obtain maximum refinement level
   maxamr = getmaxrefinement(cellid, ncells)

   vars = [node["name"] for node in n.var]

   hasvdf = let
      n.cellwithVDF[1]["arraysize"] != "0"
   end

   # File IOstream is not closed for sake of data processing later.

   meta = MetaVLSV(splitdir(file)..., fid, n, vars, celldict, cellindex,
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

function getmaxrefinement(cellid::AbstractArray{UInt, 1}, ncells::NTuple{3, UInt})
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
      for node in meta.nodeVLSV.var
         if node["name"] == var
            haskey(node, "unit") || break
            unit = node["unit"]
            unitLaTeX = node["unitLaTeX"]
            variableLaTeX = node["variableLaTeX"]
            unitConversion = node["unitConversion"]
         end
      end
   end

   VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

@inline function readcoords(fid::IOStream ,footer::EzXML.Node, qstring::String)
   node = findfirst(qstring, footer)

   arraysize = parse(Int, node["arraysize"])
   offset = parse(Int, nodecontent(node))

   # Warning: it may be Float32 in Vlasiator
   coord = Vector{Float64}(undef, arraysize)
   seek(fid, offset)
   read!(fid, coord)

   coord
end

function readvcoords(fid::IOStream, footer::EzXML.Node, species::String, qstring::String)
   local coord
   nodes = findall(qstring, footer)

   for i in eachindex(nodes)[3:end]
      if nodes[i]["mesh"] == species
         arraysize = parse(Int, nodes[i]["arraysize"])
         offset = parse(Int, nodecontent(nodes[i]))
         # Warning: it may be Float32 in Vlasiator
         coord = Vector{Float64}(undef, arraysize)
         seek(fid, offset)
         read!(fid, coord)
         break
      end
   end

   coord
end

"Return spatial mesh information."
function readmesh(fid::IOStream, footer::EzXML.Node)
   # Assume SpatialGrid and FsGrid follows Vlasiator 5 standard
   node = findfirst("//MESH_BBOX", footer)
   offset = parse(Int, nodecontent(node))

   bbox = Vector{UInt}(undef, 6)
   seek(fid, offset)
   read!(fid, bbox)

   nodeX = readcoords(fid, footer, "//MESH_NODE_CRDS_X")
   nodeY = readcoords(fid, footer, "//MESH_NODE_CRDS_Y")
   nodeZ = readcoords(fid, footer, "//MESH_NODE_CRDS_Z")

   @inbounds ncells = (bbox[1], bbox[2], bbox[3])
   @inbounds block_size = (bbox[4], bbox[5], bbox[6])
   @inbounds coordmin = (nodeX[begin], nodeY[begin], nodeZ[begin])
   @inbounds coordmax = (nodeX[end], nodeY[end], nodeZ[end])

   ncells, block_size, coordmin, coordmax
end

"Return velocity mesh information."
function readvmesh(fid::IOStream, footer::EzXML.Node, species::String)
   nodes = findall("//MESH_BBOX", footer)

   bbox = Vector{UInt}(undef, 6)

   for i in eachindex(nodes)[3:end]
      if nodes[i]["mesh"] == species
         offset = parse(Int, nodecontent(nodes[i]))
         seek(fid, offset)
         read!(fid, bbox)
         break
      end
   end

   nodeX = readvcoords(fid, footer, species, "//MESH_NODE_CRDS_X")
   nodeY = readvcoords(fid, footer, species, "//MESH_NODE_CRDS_Y")
   nodeZ = readvcoords(fid, footer, species, "//MESH_NODE_CRDS_Z")

   bbox, nodeX, nodeY, nodeZ
end

"""
    readvariable(meta::MetaVLSV, var::String, sorted::Bool=true, usemmap::Bool=false) -> Array

Return variable value of `var` from the VLSV file associated with `meta`. By default for
DCCRG variables are sorted by cell ID. `usemmap` decides whether to use memory-mapped IO,
especially for large arrays.
"""
function readvariable(meta::MetaVLSV, var::String, sorted::Bool=true, usemmap::Bool=false)
   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta)
      return v::Array
   end

   T, offset, asize, dsize, vsize = getvarinfo(meta.nodeVLSV.var, var)
   if vsize == 1
      raw = readvector(meta.fid, T, offset, asize, dsize, usemmap)
   else
      raw = readarray(meta.fid, T, offset, asize, dsize, vsize, usemmap)
   end

   if startswith(var, "fg_") # fsgrid
      bbox = @. meta.ncells * 2^meta.maxamr
      # Determine fsgrid domain decomposition
      nIORanks = readparameter(meta, "numWritingRanks")::Int32

      dataOrdered =
         if ndims(raw) > 1
            @inbounds Array{Float32}(undef, size(raw,1), bbox[1], bbox[2], bbox[3])
         else
            @inbounds Array{Float32}(undef, bbox[1], bbox[2], bbox[3])
         end

      @inbounds fgDecomposition = @views getDomainDecomposition(bbox[1:3], nIORanks)

      _fillFG!(dataOrdered, raw, fgDecomposition, nIORanks, bbox)

      v = dropdims(dataOrdered, dims=(findall(size(dataOrdered) .== 1)...,))
   elseif sorted # dccrg grid
      @inbounds v = ndims(raw) == 1 ? raw[meta.cellindex] : raw[:,meta.cellindex]
      if eltype(v) == Float64; v = Float32.(v); end
   else # dccrg grid
      v = raw
   end

   return v::Union{Array, Base.ReinterpretArray}
end

"""
    readvariable(meta::MetaVLSV, var::String, ids::Vector{<:Integer}) -> Array

Read variable `var` in a collection of cells `ids` associated with `meta`. if `ids` is
empty, return the whole sorted array of `var`.
"""
function readvariable(meta::MetaVLSV, var::String, ids::Vector{<:Integer})::Array
   startswith(var, "fg_") && error("Currently does not support reading fsgrid!")

   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta, ids)
      return v
   end

   if (local nid = length(ids)) == 0
      v = readvariable(meta, var)
   else
      T, offset, asize, dsize, vsize = getvarinfo(meta.nodeVLSV.var, var)

      v = _readcells(T, meta.fid, meta.celldict, ids, nid, offset, asize, dsize, vsize)

      if T === Float64
         v = Float32.(v)
      end
   end

   return v
end

"""
    readvariable(meta::MetaVLSV, var::String, cid::Integer) -> Array

Read variable `var` in cell `cid` associated with `meta`.
"""
function readvariable(meta::MetaVLSV, var::String, cid::Integer)::Array
   startswith(var, "fg_") && error("Currently does not support reading fsgrid!")

   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta, cid)
      return v
   end

   T, offset, _, dsize, vsize = getvarinfo(meta.nodeVLSV.var, var)

   v = _readcell(T, meta.fid, meta.celldict, cid, offset, dsize, vsize)

   if T === Float64
      v = Float32.(v)
   end

   return v
end

@inline function _readcells(::Type{T}, fid::IOStream,
   celldict::Dict{UInt, Int}, ids::Vector{<:Integer},
   nid::Int, offset::Int, asize::Int, dsize::Int, vsize::Int)::Array{T} where T
   v = vsize == 1 ? Vector{T}(undef, nid) : Array{T,2}(undef, vsize, nid)

   w = let
      a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
      reshape(reinterpret(T, a), vsize, asize)
   end

   _fillv!(v, w, celldict, ids)

   v
end

@inline function _readcell(::Type{T}, fid::IOStream, celldict::Dict{UInt, Int},
   cid::Integer, offset::Int, dsize::Int, vsize::Int)::Array{T} where T

   v = vsize == 1 ? Vector{T}(undef, 1) : Array{T,2}(undef, vsize, 1)
   seek(fid, offset + (celldict[cid]-1)*vsize*dsize)
   read!(fid, v)
end

@inline function _fillv!(v::Vector{T}, w::AbstractArray{T}, celldict::Dict{UInt, Int},
   ids::Vector{<:Integer}) where T

   for i in eachindex(ids)
      @inbounds v[i] = w[celldict[ids[i]]]
   end
end

@inline function _fillv!(v::Matrix{T}, w::AbstractArray{T}, celldict::Dict{UInt, Int},
   ids::Vector{<:Integer}) where T

   for i in eachindex(ids), iv in axes(v,1)
      @inbounds v[iv,i] = w[iv,celldict[ids[i]]]
   end
end

function _fillFG!(dataOrdered::Array{Float32}, raw::AbstractArray{<:Real},
   fgDecomposition::NTuple{3, Int}, nIORanks::Int32, bbox::NTuple{3, Int})
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

@inline function getcellid(fid::IOStream, nodevar::AbstractVector{EzXML.Node})
   _, offset, asize, _, _ = getvarinfo(nodevar, "CellID")
   cellid = mmap(fid, Vector{UInt}, asize, offset)
end

"""
    extractsat(files::AbstractVector{String}, var::String, cid::Integer)

Extract `var` at a fixed cell ID `cid` from `files`. This assumes that `files` come from the
same grid structure.
"""
function extractsat(files::AbstractVector{String}, var::String, cid::Integer)
   v = open(files[1], "r") do fid
      footer = getfooter(fid)
      nodevar = findall("//VARIABLE", footer)
      T, _, _, _, vsize = getvarinfo(nodevar, var)

      Array{T,2}(undef, vsize, length(files))
   end

   for i in eachindex(files)
      fid = open(files[i], "r")
      footer = getfooter(fid)
      nodevar = findall("//VARIABLE", footer)
      cellid = getcellid(fid, nodevar)
      c_ = findfirst(isequal(cid), cellid)
      _, offset, _, dsize, vsize = getvarinfo(nodevar, var)
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
    readparameter(meta::MetaVLSV, param::String)

Return the parameter value from the VLSV file associated with `meta`.
"""
readparameter(meta::MetaVLSV, param::String) =
   readparameter(meta.fid, meta.nodeVLSV.param, param)

function readparameter(fid::IOStream, nodeparam::NodeVector, name::String)
   T, offset = getparaminfo(nodeparam, name)
   seek(fid, offset)
   p = read(fid, T)
end

"""
    hasparameter(meta::MetaVLSV, param::String) -> Bool

Check if the VLSV file contains a certain parameter `param`.
"""
hasparameter(meta::MetaVLSV, param::String) = hasname(meta.nodeVLSV.param, param)

"""
    hasvariable(meta::MetaVLSV, var::String) -> Bool

Check if the VLSV file associated with `meta` contains a variable `var`.
"""
hasvariable(meta::MetaVLSV, var::String) = hasname(meta.nodeVLSV.var, var)

"Check if the XML `nodes` contain a node with `name`."
function hasname(nodes::NodeVector, name::String)
   any(node -> node["name"] == name, nodes)
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
   (;fid, nodeVLSV) = meta
   (;vblock_size) = meta.meshes[species]
   bsize = prod(vblock_size)

   local offset_v::Int, nblocks::Int

   let cellsWithVDF, nblock_C
      for node in nodeVLSV.cellwithVDF
         if node["name"] == species
            asize = parse(Int, node["arraysize"])
            offset = parse(Int, nodecontent(node))
            cellsWithVDF = Vector{UInt}(undef, asize)
            seek(fid, offset)
            read!(fid, cellsWithVDF)
            break
         end
      end

      for node in nodeVLSV.cellblocks
         if node["name"] == species
            asize = parse(Int, node["arraysize"])
            dsize = parse(Int, node["datasize"])
            offset = parse(Int, nodecontent(node))
            nblock_C = dsize == 4 ?
               Vector{UInt32}(undef, asize) : Vector{UInt}(undef, asize)
            seek(fid, offset)
            read!(fid, nblock_C)
            break
         end
      end

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
      offset_v = sum(@view nblock_C[1:cellWithVDFIndex-1])
   end

   local dsize, vsize, offset
   # Read in avgs
   for node in nodeVLSV.blockvar
      if node["name"] == species
         dsize = parse(Int, node["datasize"])
         vsize = parse(Int, node["vectorsize"])
         offset = parse(Int, nodecontent(node))
         break
      end
   end

   data = let
      Tavg = dsize == 4 ? Float32 : Float64
      a = mmap(fid, Vector{UInt8}, dsize*vsize*nblocks, offset_v*vsize*dsize + offset)
      reshape(reinterpret(Tavg, a), vsize, nblocks)
   end

   # Read in block IDs
   for node in nodeVLSV.blockid
      if node["name"] == species
         dsize = parse(Int, node["datasize"])
         offset = parse(Int, nodecontent(node))
         break
      end
   end

   blockIDs = let T = dsize == 4 ? UInt32 : UInt
      Vector{T}(undef, nblocks)
   end
   seek(fid, offset_v*dsize + offset)
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