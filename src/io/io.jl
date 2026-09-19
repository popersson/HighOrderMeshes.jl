###########################################################################
## Binary .hom format
#
# A .hom file is a sequence of named arrays ("records") stored as raw
# little-endian binary. It is self-describing and readable from any language
# with a few lines of code:
#
#   [8 bytes]         magic string "HOMESHv1"
#   [Int64]           number of records
#   then, for each record:
#     [Int64]         length of the name in bytes
#     [bytes]         name (ASCII)
#     [Int64]         element type code (see hom_type_codes)
#     [Int64]         number of dimensions (0 for a scalar)
#     [Int64 × ndims] shape
#     [data]          elements in column-major order
#
# Type codes: 1 Float64, 2 Float32, 3 Float16, 11 Int64, 12 Int32, 13 Int16, 14 UInt8.
#
# Records written by savemesh (D = spatial dimension, ns = nodes per element,
# nf = faces per element, T = coordinate type):
#
#   elgeom     Int64  scalar        1 = Simplex, 2 = Block
#   ref_nodes  T      ns × D        reference nodes of the FiniteElement
#   x          T      nnodes × D    node coordinates
#   el         Int64  ns × nel      element-to-node table (one-based)
#   nb_el      Int32  nf × nel      neighbor element (<= 0: boundary face with tag -value)
#   nb_face    Int16  nf × nel      neighbor face index
#   nb_perm    Int16  nf × nel      neighbor face permutation
#
# The polynomial degree follows from the number of reference nodes. Readers
# must ignore records they do not know, so that records can be added later.

const hom_magic = "HOMESHv1"
const hom_type_codes = Dict{Int64,DataType}(
    1 => Float64, 2 => Float32, 3 => Float16,
    11 => Int64, 12 => Int32, 13 => Int16, 14 => UInt8)
const hom_code_of = Dict{DataType,Int64}(v => k for (k, v) in hom_type_codes)

const _little_endian = Base.ENDIAN_BOM == 0x04030201

function write_record(io::IO, name::AbstractString, a::AbstractArray)
    code = get(hom_code_of, eltype(a), nothing)
    isnothing(code) && error("Unsupported element type $(eltype(a)) for .hom record \"$name\"")
    write(io, htol(Int64(ncodeunits(name))))
    write(io, name)
    write(io, htol(Int64(code)))
    write(io, htol(Int64(ndims(a))))
    for s in size(a)
        write(io, htol(Int64(s)))
    end
    data = a isa Array ? a : collect(a)
    write(io, _little_endian ? data : htol.(data))
    nothing
end

function read_record(io::IO)
    n    = ltoh(read(io, Int64))
    name = String(read(io, n))
    code = ltoh(read(io, Int64))
    T    = get(hom_type_codes, code, nothing)
    isnothing(T) && error("Unknown element type code $code for .hom record \"$name\"")
    nd    = Int(ltoh(read(io, Int64)))
    shape = ntuple(_ -> Int(ltoh(read(io, Int64))), nd)
    a = Array{T,nd}(undef, shape)
    read!(io, a)
    _little_endian || (a .= ltoh.(a))
    name, a
end

"""
    savemesh(fname, m::HighOrderMesh)

Save `m` to a binary `.hom` file. The format stores everything needed to
reconstruct the mesh exactly: geometry type, reference nodes, coordinates,
connectivity and neighbor data, as named little-endian arrays (see the
header of `io.jl` for the layout). Use `loadmesh` to read it back.
"""
function savemesh(fname, m::HighOrderMesh{D,G,T}) where {D,G,T}
    geo_id = findfirst(==(Base.typename(G).wrapper), geometry_types)
    isnothing(geo_id) && error("Unsupported geometry type: $G. Add to geometry_types.")
    records = [
        ("elgeom",    fill(Int64(geo_id))),
        ("ref_nodes", ref_nodes(m.fe)),
        ("x",         m.x),
        ("el",        Matrix{Int64}(m.el)),
        ("nb_el",     getindex.(m.nb, 1)),
        ("nb_face",   getindex.(m.nb, 2)),
        ("nb_perm",   getindex.(m.nb, 3)),
    ]
    open(fname, "w") do io
        write(io, hom_magic)
        write(io, htol(Int64(length(records))))
        for (name, a) in records
            write_record(io, name, a)
        end
    end
    nothing
end

"""
    loadmesh(fname) -> HighOrderMesh

Load a `HighOrderMesh` from a binary `.hom` file previously written by `savemesh`.
"""
function loadmesh(fname)
    open(fname, "r") do io
        magic = String(read(io, ncodeunits(hom_magic)))
        magic == hom_magic || error("Not a .hom file: bad magic string \"$magic\"")
        nrec = ltoh(read(io, Int64))
        rec  = Dict{String,Any}()
        for _ in 1:nrec
            name, a = read_record(io)
            rec[name] = a
        end
        for key in ("elgeom", "ref_nodes", "x", "el", "nb_el", "nb_face", "nb_perm")
            haskey(rec, key) || error("Missing record \"$key\" in .hom file")
        end

        x  = rec["x"]
        D  = size(x, 2)
        T  = eltype(x)
        eg = geometry_types[rec["elgeom"][]]{D}()
        fe = FiniteElement(eg, Matrix{T}(rec["ref_nodes"]))
        el = Matrix{Int}(rec["el"])
        nb = tuple.(Int32.(rec["nb_el"]), Int16.(rec["nb_face"]), Int16.(rec["nb_perm"]))
        HighOrderMesh{D,typeof(eg),T}(fe, x, el, nb)
    end
end

###########################################################################
## Text helpers (used by the VTK writer)

# Write matrix `x` to an open IO stream, one row per line, space-separated.
function write_matrix(f, x)
    for i = 1:size(x,1)
        join(f, x[i,:], ' ')
        println(f)
    end
end
