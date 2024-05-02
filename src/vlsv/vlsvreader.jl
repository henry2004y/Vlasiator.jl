# VLSV reader in Julia

const NodeVector = SubArray{Node, 1, Vector{Node}, Tuple{UnitRange{Int64}},
   true}

@lazy struct VMeshInfo
   "number of velocity blocks"
   vblocks::NTuple{3, Int}
   vblocksize::NTuple{3, Int}
   vmin::NTuple{3, Float64}
   vmax::NTuple{3, Float64}
   dv::NTuple{3, Float64}
   @lazy cellwithVDF::Vector{Int}
   @lazy nblock_C::Vector{Int}
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
   celldict::Dict{Int, Int}
   "ordered sequence indexes of raw cell IDs"
   cellindex::Vector{Int}
   time::Float64
   maxamr::Int
   hasvdf::Bool
   ncells::NTuple{3, Int}
   blocksize::NTuple{3, Int}
   coordmin::NTuple{3, Float64}
   coordmax::NTuple{3, Float64}
   dcoord::NTuple{3, Float64}
   species::Vector{String}
   meshes::Dict{String, VMeshInfo}
end


function Base.show(io::IO, meta::MetaVLSV)
   print(io, "File       : ")
   printstyled(io, meta.name, '\n'; color=:cyan, underline=true)
   print(io, "Time       : ")
   t = @sprintf "%.4e s" meta.time
   printstyled(io, t, '\n'; color=:cyan)
   print(io, "Dimension  : ")
   printstyled(io, ndims(meta), '\n'; color=:yellow)
   print(io, "Max AMR lvl: ")
   printstyled(io, meta.maxamr, '\n'; color=:yellow)
   print(io, "Has VDF    : ")
   printstyled(io, meta.hasvdf, '\n'; color=:yellow)
   print(io, "Variables  : ")
   printstyled(io, meta.variable, '\n'; color=:green)
end

function Base.show(io::IO, s::VarInfo)
   println(io, "Variable in LaTeX: ", s.variableLaTeX)
   println(io, "Unit: ", s.unit)
   println(io, "Unit in LaTeX: ", s.unitLaTeX)
   println(io, "Unit conversion: ", s.unitConversion)
end

function Base.show(io::IO, vmesh::VMeshInfo)
   println(io, "vblocks: ", vmesh.vblocks)
   println(io, "vblock size: ", vmesh.vblocksize)
   foreach((vmin,dv,vmax,comp) -> println(io, "v$comp range: ", vmin, ":", dv, ":", vmax),
      vmesh.vmin, vmesh.dv, vmesh.vmax, 'x':'z')
end

"Return the XML footer of opened VLSV file."
@inline function getfooter(fid::IOStream)
   # First 8 bytes indicate big-endian or else
   endian_offset = 8
   seek(fid, endian_offset)
   # Obtain the offset of the XML
   offset = read(fid, Int)
   seek(fid, offset)
   str = read(fid, String)
   footer = parse(Node, str)
end

@inline function getdatatype(datatype::Symbol, datasize::Int)
   T::DataType =
      if datatype == :float
         datasize == 4 ? Float32 : Float64
      elseif datatype == :int
         datasize == 4 ? Int32 : Int
      elseif datatype == :uint
         datasize == 4 ? UInt32 : UInt
      else
         throw(ArgumentError("unknown type $datatype"))
      end
end

function getvarinfo(nodevar::NodeVector, name::String)
   local arraysize, datasize, datatype, vectorsize, offset
   isFound = false

   for nv in nodevar
      var = attributes(nv)
      if var["name"] == name
         arraysize = Parsers.parse(Int, var["arraysize"])
         datasize = Parsers.parse(Int, var["datasize"])
         datatype = Symbol(var["datatype"])
         vectorsize = Parsers.parse(Int, var["vectorsize"])
         offset = Parsers.parse(Int, value(nv[1]))
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
      param = attributes(p)
      if param["name"] == name
         datasize = Parsers.parse(Int, param["datasize"])
         datatype = Symbol(param["datatype"])
         offset = Parsers.parse(Int, value(p[1]))
         isFound = true
         break
      end
   end

   isFound || throw(ArgumentError("unknown variable $name"))

   T = getdatatype(datatype, datasize)

   T, offset
end

function Base.findall(xpath::AbstractString, nodes::Vector{Node})
   findall(x -> tag(x) == xpath, nodes)
end

function Base.findfirst(xpath::AbstractString, nodes::Vector{Node})
   findfirst(x -> tag(x) == xpath, nodes)
end

"General inquiry of element `tag` with `tagname` and `attr`."
function getObjInfo(ns::Vector{Node}, name::String, tagname::String, attr::String)
   local arraysize, datasize, datatype, vectorsize, offset
   isFound = false

   for i in findall(tagname, ns)
      var = attributes(ns[i])
      if var[attr] == name
         arraysize = Parsers.parse(Int, var["arraysize"])
         datasize = Parsers.parse(Int, var["datasize"])
         datatype = Symbol(var["datatype"])
         vectorsize = Parsers.parse(Int, var["vectorsize"])
         offset = Parsers.parse(Int, value(ns[i][1]))
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
   ::Val{false}) where T
   w = Vector{T}(undef, asize)
   seek(fid, offset)
   read!(fid, w)

   w
end

@inline function readvector(fid::IOStream, ::Type{T}, offset::Int, asize::Int, dsize::Int,
   ::Val{true}) where T
   w =
      if offset % dsize == 0
         mmap(fid, Vector{T}, asize, offset)
      else
         a = mmap(fid, Vector{UInt8}, dsize*asize, offset)
         reinterpret(T, a)
      end
end

@inline function readarray(fid::IOStream, ::Type{T}, offset::Int, asize::Int, dsize::Int,
   vsize::Int, ::Val{false}) where T
   w = Array{T,2}(undef, vsize, asize)
   seek(fid, offset)
   read!(fid, w)

   w
end

@inline function readarray(fid::IOStream, ::Type{T}, offset::Int, asize::Int, dsize::Int,
   vsize::Int, ::Val{true}) where T
   w =
      if offset % dsize == 0
         mmap(fid, Array{T,2}, (vsize, asize), offset)
      else
         a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
         reshape(reinterpret(T, a), vsize, asize)
      end
end

"""
    load(file::AbstractString)) -> MetaVLSV

Generate a MetaVLSV object from `file` of VLSV format.
"""
function load(file::AbstractString)
   fid = open(file, "r")

   footer = getfooter(fid)
   ns = children(footer[end])

   @no_escape begin
      ibegin_, iend_ = @alloc(Int, 6), @alloc(Int, 6)
      for i in eachindex(ibegin_)
         ibegin_[i] = 0
      end

      for i in eachindex(ns)
         if tag(ns[i]) == "VARIABLE"
            if ibegin_[1] == 0 ibegin_[1] = i end
            iend_[1] = i
         elseif tag(ns[i]) == "PARAMETER"
            if ibegin_[2] == 0 ibegin_[2] = i end
            iend_[2] = i
         elseif tag(ns[i]) == "CELLSWITHBLOCKS"
            if ibegin_[3] == 0 ibegin_[3] = i end
            iend_[3] = i
         elseif tag(ns[i]) == "BLOCKSPERCELL"
            if ibegin_[4] == 0 ibegin_[4] = i end
            iend_[4] = i
         elseif tag(ns[i]) == "BLOCKVARIABLE"
            if ibegin_[5] == 0 ibegin_[5] = i end
            iend_[5] = i
         elseif tag(ns[i]) == "BLOCKIDS"
            if ibegin_[6] == 0 ibegin_[6] = i end
            iend_[6] = i
         end
      end

      @views begin
         nodevar = ns[ibegin_[1]:iend_[1]]
         nodeparam = ns[ibegin_[2]:iend_[2]]
         if ibegin_[3] != 0
            nodecellwithVDF = ns[ibegin_[3]:iend_[3]]
            nodecellblocks = ns[ibegin_[4]:iend_[4]]
            nodeblockvar = ns[ibegin_[5]:iend_[5]]
            nodeblockid = ns[ibegin_[6]:iend_[6]]
         else # non-standard Vlasiator outputs with no vspace meta data
            nodecellwithVDF = ns[1:0]
            nodecellblocks = ns[1:0]
            nodeblockvar = ns[1:0]
            nodeblockid = ns[1:0]
         end
      end
   end

   n = NodeVLSV(nodevar, nodeparam, nodecellwithVDF, nodecellblocks, nodeblockvar,
      nodeblockid)

   cellid = getcellid(fid, n.var)
   cellindex = sortperm(cellid, alg=MergeSort)

   ncells, blocksize, coordmin, coordmax = readmesh(fid, ns)

   dcoord = ntuple(i -> (coordmax[i] - coordmin[i]) / ncells[i], Val(3))

   meshes = Dict{String, VMeshInfo}()

   # Find all species by the BLOCKIDS tag
   species = String[]

   @no_escape begin
      vblocks = @alloc(Int, 3)
      vblocksize = @alloc(Int, 3)
      vmin = @alloc(Float64, 3)
      vmax = @alloc(Float64, 3)

      for node in n.blockid
         at = attributes(node)
         if haskey(at, "name")
            # VLSV 5.0 file with bounding box
            popname = at["name"]

            vbox, nodeX, nodeY, nodeZ = readvmesh(fid, ns, popname)

            vblocks[:] = @view vbox[1:3]
            vblocksize[:] = @view vbox[4:6]
            vmin[1], vmin[2], vmin[3] = nodeX[begin], nodeY[begin], nodeZ[begin]
            vmax[1], vmax[2], vmax[3] = nodeX[end], nodeY[end], nodeZ[end]
         else
            # In VLSV before 5.0 the mesh is defined with parameters.
            popname = "avgs"
            if "vxblocks_ini" in getindex.(n.param, "name")
               vblocks_str = ("vxblocks_ini", "vyblocks_ini", "vzblocks_ini")
               vmin_str = ("vxmin", "vymin", "vzmin")
               vmax_str = ("vxmax", "vymax", "vzmax")
               map!(i -> readparameter(fid, n.param, vblocks_str[i]), vblocks, 1:3)
               map!(i -> readparameter(fid, n.param, vmin_str[i]), vmin, 1:3)
               map!(i -> readparameter(fid, n.param, vmax_str[i]), vmax, 1:3)
               vblocksize .= 4
            else
               error("File not written by Vlasiator!")
            end
         end

         dv = ntuple(i -> (vmax[i] - vmin[i]) / vblocks[i] / vblocksize[i], Val(3))

         # Update list of active species
         if popname ∉ species
            push!(species, popname)
         end

         # Create a new object for this species
         popVMesh = VMeshInfo(NTuple{3}(vblocks), NTuple{3}(vblocksize),
            NTuple{3}(vmin), NTuple{3}(vmax), dv, uninit, uninit)

         meshes[popname] = popVMesh
      end
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

   vars = [attributes(node)["name"] for node in n.var]

   hasvdf = let
      if length(n.cellwithVDF) == 0
         false
      else
         attributes(n.cellwithVDF[1])["arraysize"] != "0"
      end
   end

   # File IOstream is not closed for sake of data processing later.

   meta = MetaVLSV(splitdir(file)..., fid, n, vars, celldict, cellindex,
      timesim, maxamr, hasvdf, ncells, blocksize, coordmin, coordmax,
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

"Get maximum refinement level, assuming 3D grid and dichotomy method."
function getmaxrefinement(cellid::Vector{Int}, ncells::NTuple{3, Int})
   ncell = prod(ncells)
   maxamr, cid = 0, ncell
   while cid < maximum(cellid)
      maxamr += 1
      cid += ncell << (3*maxamr)
   end

   maxamr
end

"""
    readvariablemeta(meta, var) -> VarInfo

Return VarInfo of `var` in the VLSV file associated with `meta`.
"""
function readvariablemeta(meta::MetaVLSV, var::String)
   varSym = isa(var, AbstractString) ? Symbol(var) : var

   unit, unitLaTeX, variableLaTeX, unitConversion = "", "", "", ""

   if varSym in keys(units_predefined)
      unit, variableLaTeX, unitLaTeX = units_predefined[varSym]
   elseif hasvariable(meta, var) # For Vlasiator 5 files, MetaVLSV is included
      for node in meta.nodeVLSV.var
         at = attributes(node)
         if at["name"] == var && haskey(at, "unit")
            unit = at["unit"]
            unitLaTeX = at["unitLaTeX"]
            variableLaTeX = at["variableLaTeX"]
            unitConversion = at["unitConversion"]
            break
         end
      end
   end

   VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end

@inline function readcoords(fid::IOStream, ns::Vector{Node}, qstring::String)
   node = ns[findfirst(qstring, ns)]

   arraysize = Parsers.parse(Int, attributes(node)["arraysize"])
   offset = Parsers.parse(Int, value(node[1]))

   # Warning: it may be Float32 in Vlasiator
   coord = Vector{Float64}(undef, arraysize)
   seek(fid, offset)
   read!(fid, coord)

   coord
end

function readvcoords(fid::IOStream, ns::Vector{Node}, species::String, qstring::String)
   local coord
   is = findall(qstring, ns)

   for i in reverse(is)
      at = attributes(ns[i])
      if at["mesh"] == species
         arraysize = Parsers.parse(Int, at["arraysize"])
         offset = Parsers.parse(Int, value(ns[i][1]))
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
function readmesh(fid::IOStream, ns::Vector{Node})
   # Assume SpatialGrid and FsGrid follows Vlasiator 5 standard
   node = ns[findfirst("MESH_BBOX", ns)]
   offset = Parsers.parse(Int, value(node[1]))

   bbox = Vector{Int}(undef, 6)
   seek(fid, offset)
   read!(fid, bbox)

   nodeX = readcoords(fid, ns, "MESH_NODE_CRDS_X")
   nodeY = readcoords(fid, ns, "MESH_NODE_CRDS_Y")
   nodeZ = readcoords(fid, ns, "MESH_NODE_CRDS_Z")

   @inbounds ncells = (bbox[1], bbox[2], bbox[3])
   @inbounds blocksize = (bbox[4], bbox[5], bbox[6])
   @inbounds coordmin = (nodeX[begin], nodeY[begin], nodeZ[begin])
   @inbounds coordmax = (nodeX[end], nodeY[end], nodeZ[end])

   ncells, blocksize, coordmin, coordmax
end

"Return velocity mesh information."
function readvmesh(fid::IOStream, ns::Vector{Node}, species::String)
   is = findall("MESH_BBOX", ns)

   bbox = Vector{Int}(undef, 6)

   for i in reverse(is)
      at = attributes(ns[i])
      if at["mesh"] == species
         offset = Parsers.parse(Int, value(ns[i][1]))
         seek(fid, offset)
         read!(fid, bbox)
         break
      end
   end

   nodeX = readvcoords(fid, ns, species, "MESH_NODE_CRDS_X")
   nodeY = readvcoords(fid, ns, species, "MESH_NODE_CRDS_Y")
   nodeZ = readvcoords(fid, ns, species, "MESH_NODE_CRDS_Z")

   bbox, nodeX, nodeY, nodeZ
end

"""
    readvariable(meta::MetaVLSV, var::String, sorted::Bool=true, usemmap::Bool=false) -> Array

Return variable value of `var` from the VLSV file associated with `meta`. By default for
DCCRG variables are sorted by cell ID. `usemmap` decides whether to use memory-mapped IO,
which is especially useful for large arrays.
"""
function readvariable(meta::MetaVLSV, var::String, sorted::Bool=true, usemmap::Bool=false)
   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta)
      return v::Array
   end

   T, offset, asize, dsize, vsize = getvarinfo(meta.nodeVLSV.var, var)
   if vsize == 1
      raw = readvector(meta.fid, T, offset, asize, dsize, Val(usemmap))
   else
      raw = readarray(meta.fid, T, offset, asize, dsize, vsize, Val(usemmap))
   end

   if startswith(var, "fg_") # fsgrid
      bbox = @. meta.ncells * 2^meta.maxamr
      # Determine fsgrid domain decomposition
      nIORanks = readparameter(meta, "numWritingRanks")::Int32

      v =
         if ndims(raw)::Int64 > 1
            @inbounds Array{Float32}(undef, size(raw,1)::Int64, bbox...)
         else
            @inbounds Array{Float32}(undef, bbox...)
         end

      fgDecomposition = getDomainDecomposition(bbox, nIORanks)

      _fillFG!(v, raw, fgDecomposition, nIORanks, bbox)
   elseif sorted # dccrg grid
      @inbounds v = ndims(raw) == 1 ? raw[meta.cellindex] : raw[:,meta.cellindex]
      if T === Float64; v = map(x->Float32(x), v); end
   else # dccrg grid
      v = raw
   end

   return v::Union{Array, Base.ReinterpretArray}
end

"""
    readvariable(meta::MetaVLSV, var::String, ids::Vector{Int}) -> Array

Read variable `var` in a collection of cells `ids` associated with `meta`. if `ids` is
empty, return the whole sorted array of `var`.
"""
function readvariable(meta::MetaVLSV, var::String, ids::Vector{Int})::Array
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
         v = map(x->Float32(x), v)
      end
   end

   return v
end

"""
    readvariable(meta::MetaVLSV, var::String, cid::Int) -> Array

Read variable `var` in cell `cid` associated with `meta`.
"""
function readvariable(meta::MetaVLSV, var::String, cid::Int)::Array
   startswith(var, "fg_") && error("Currently does not support reading fsgrid!")

   if (local symvar = Symbol(var)) in keys(variables_predefined)
      v = variables_predefined[symvar](meta, [cid])
      return v
   end

   T, offset, _, dsize, vsize = getvarinfo(meta.nodeVLSV.var, var)

   v = _readcell(T, meta.fid, meta.celldict, cid, offset, dsize, vsize)

   if T === Float64
      v = map(x->Float32(x), v)
   end

   return v
end

@inline function _readcells(::Type{T}, fid::IOStream,
   celldict::Dict{Int, Int}, ids::Vector{Int},
   nid::Int, offset::Int, asize::Int, dsize::Int, vsize::Int)::Array{T} where T
   v = vsize == 1 ? Vector{T}(undef, nid) : Array{T,2}(undef, vsize, nid)

   if offset % dsize == 0
      w = mmap(fid, Array{T,2}, (vsize, asize), offset)
   else
      a = mmap(fid, Vector{UInt8}, dsize*vsize*asize, offset)
      w = reshape(reinterpret(T, a), vsize, asize)
   end

   _fillv!(v, w, celldict, ids)

   v
end

@inline function _readcell(::Type{T}, fid::IOStream, celldict::Dict{Int, Int},
   cid::Int, offset::Int, dsize::Int, vsize::Int)::Array{T} where T

   v = vsize == 1 ? Vector{T}(undef, 1) : Array{T,2}(undef, vsize, 1)
   seek(fid, offset + (celldict[cid]-1)*vsize*dsize)
   read!(fid, v)
end

@inline function _fillv!(v::Vector{T}, w::AbstractArray{T}, celldict::Dict{Int, Int},
   ids::Vector{Int}) where T

   for i in eachindex(ids)
      @inbounds v[i] = w[celldict[ids[i]]]
   end
end

@inline function _fillv!(v::Matrix{T}, w::AbstractArray{T}, celldict::Dict{Int, Int},
   ids::Vector{Int}) where T

   for i in eachindex(ids), iv in axes(v,1)
      @inbounds v[iv,i] = w[iv,celldict[ids[i]]]
   end
end

function _fillFG!(dataOrdered::Array{Float32}, raw::AbstractArray{<:Real},
   fgDecomposition::NTuple{3, Int}, nIORanks::Int32, bbox::NTuple{3, Int})
   offsetnow = 1

   @inbounds @views for i in 1:nIORanks
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

@inline function getcellid(fid::IOStream, nodevar::NodeVector)
   _, offset, asize, _, _ = getvarinfo(nodevar, "CellID")
   cellid = mmap(fid, Vector{Int}, asize, offset)
end

"""
    extractsat(files::AbstractVector{String}, var::String, cid::Int)

Multi-threaded extraction of `var` at a fixed cell ID `cid` from `files`. This assumes that
`files` come from the same grid structure.
"""
function extractsat(files::AbstractVector{String}, var::String, cid::Int)
   v = open(files[1], "r") do fid
      footer = getfooter(fid)
      ns = children(footer[end])
      is = findall("VARIABLE", ns)
      nodevar = @view ns[is[1]:is[end]]
      T, _, _, _, vsize = getvarinfo(nodevar, var)

      Array{T,2}(undef, vsize, length(files))
   end

   Threads.@threads for i in eachindex(files)
      fid = open(files[i], "r")
      footer = getfooter(fid)
      ns = children(footer[end])
      is = findall("VARIABLE", ns)
      nodevar = @view ns[is[1]:is[end]]
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

   data = ndims(data) == 1 || size(data, 1) == 1 || ndims(data) == 3 ?
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

   @inbounds for i in 1:min(nprocs, globalsize[1])
      iBox = max(globalsize[1]/i, 1.0)

      for j in 1:min(nprocs, globalsize[2])
         i * j > nprocs && break

         jBox = max(globalsize[2]/j, 1.0)

         for k in 1:min(nprocs, globalsize[3])
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

"Check if the XML nodes `ns` contain a node of `name`."
hasname(ns::NodeVector, name::String) = any(n -> attributes(n)["name"] == name, ns)

"""
    ndims(meta::MetaVLSV) -> Int

Return the simulation spatial dimension of VLSV data.
"""
Base.ndims(meta::MetaVLSV) = count(>(1), meta.ncells)

"Determine whether `meta` is not yet closed."
Base.isopen(meta::MetaVLSV) = isopen(meta.fid)

"""
    readvcells(meta::MetaVLSV, cid::Int; species="proton") -> vcellids, vcellf

Read velocity cells of `species` from a spatial cell of ID `cid` associated with `meta`, and
return a map of velocity cell ids `vcellids` and corresponding value `vcellf`.
"""
function readvcells(meta::MetaVLSV, cid::Int; species::String="proton")
   (;fid, nodeVLSV) = meta
   (;vblocksize) = meta.meshes[species]

   init_cellswithVDF!(meta, species)

   cellWithVDFIndex = findfirst(==(cid), meta.meshes[species].cellwithVDF)
   if isnothing(cellWithVDFIndex)
      throw(ArgumentError("Cell ID $cid does not store velocity distribution!"))
   end
   @inbounds nblocks = meta.meshes[species].nblock_C[cellWithVDFIndex]
   # Offset position to vcell storage
   offset_v = sum(@view meta.meshes[species].nblock_C[1:cellWithVDFIndex-1]; init=0)

   local dsize, offset
   # Read raw VDF
   for node in nodeVLSV.blockvar
      at = attributes(node)
      if at["name"] == species
         dsize = Parsers.parse(Int, at["datasize"])
         offset = Parsers.parse(Int, value(node[1]))
         break
      end
   end

   bsize = prod(vblocksize)

   data = let
      T = dsize == 4 ? Float32 : Float64
      offset += offset_v*bsize*dsize
      a =
         if offset % dsize == 0
            mmap(fid, Array{T, 2}, (bsize, nblocks), offset)
         else
            aflat = mmap(fid, Vector{UInt8}, dsize*bsize*nblocks, offset)
            reshape(reinterpret(T, aflat), bsize, nblocks)
         end
   end

   # Read block IDs
   for node in nodeVLSV.blockid
      at = attributes(node)
      if at["name"] == species
         dsize = Parsers.parse(Int, at["datasize"])
         offset = Parsers.parse(Int, value(node[1]))
         break
      end
   end

   blockIDs = let T = dsize == 4 ? Int32 : Int
      Vector{T}(undef, nblocks)
   end
   seek(fid, offset_v*dsize + offset)
   read!(fid, blockIDs)

   # Velocity cell IDs and distributions (ordered by blocks)
   vcellids = Vector{Int32}(undef, bsize*nblocks)
   vcellf = Vector{Float32}(undef, bsize*nblocks)

   _fillvcell!(vcellids, vcellf, data, blockIDs, bsize)

   vcellids, vcellf
end

"Lazily initialize arrays of cell ID with VDF and number of vblock per cell."
function init_cellswithVDF!(meta::MetaVLSV, species::String)
   (; fid, nodeVLSV) = meta
   mesh = meta.meshes[species]

   if !@isinit mesh.cellwithVDF
      for node in nodeVLSV.cellwithVDF
         at = attributes(node)
         if at["name"] == species
            asize = Parsers.parse(Int, at["arraysize"])
            offset = Parsers.parse(Int, value(node[1]))
            cellwithVDF = Vector{Int}(undef, asize)
            seek(fid, offset)
            read!(fid, cellwithVDF)
            @init! mesh.cellwithVDF = cellwithVDF
            break
         end
      end

      for node in nodeVLSV.cellblocks
         at = attributes(node)
         if at["name"] == species
            asize = Parsers.parse(Int, at["arraysize"])
            dsize = Parsers.parse(Int, at["datasize"])
            offset = Parsers.parse(Int, value(node[1]))
            T = dsize == 4 ? Int32 : Int
            nblock_C = Vector{Int}(undef, asize)
            seek(fid, offset)
            @inbounds for i in 1:asize
               nblock_C[i] = read(fid, T)
            end
            @init! mesh.nblock_C = nblock_C
            break
         end
      end
   end

   innerBCCells = findall(==(0), mesh.nblock_C)
   deleteat!(mesh.cellwithVDF, innerBCCells)
end

@inline function _fillvcell!(vcellids, vcellf, data, blockIDs, bsize)
   @inbounds for i in eachindex(blockIDs), j in 1:bsize
      index_ = (i-1)*bsize + j
      vcellids[index_] = j - 1 + bsize*blockIDs[i]
      vcellf[index_] = data[j,i]
   end
end