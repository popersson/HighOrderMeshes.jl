###########################################################################
## Binary .hom format
#
# File layout (version "HOM_v1"):
#   [6 bytes]  magic header "HOM_v1"
#   [Int64]    D  (spatial dimension)
#   [Int64]    P  (polynomial order)
#   [UInt8]    geometry index into geometry_types
#   [UInt8]    type index into coordinate_types
#   [Int64 + data]  1D reference nodes (length + values)
#   [2×Int64 + data] x  (nnodes × D coordinate matrix)
#   [2×Int64 + data] el (nnodes_per_elem × nelems, stored as Int64)
#   [2×Int64 + data] nb (nfaces × nelems NeighborData matrix)
#
# geometry_types:   [Simplex, Block]
# coordinate_types: [Float64, Float32, Float16, BigFloat]

const coordinate_types = [Float64, Float32, Float16, BigFloat]

"""
    savemesh(fname, m::HighOrderMesh)

Save `m` to a binary `.hom` file. The format stores all mesh data needed to
exactly reconstruct the mesh, including geometry type, polynomial order,
reference nodes, coordinates, connectivity, and neighbor data.
Use `loadmesh` to read it back.
"""
function savemesh(fname, m::HighOrderMesh{D,G,P,T}) where {D,G,P,T}
    open(fname, "w") do io
        write(io, "HOM_v1") # Magic Header

        # --- 1. METADATA (D, P, Geometry) ---
        write(io, Int64(D))
        write(io, Int64(P))

        # Identify Geometry: Block{2} -> Block
        base_geo = Base.typename(G).wrapper
        geo_id = findfirst(==(base_geo), geometry_types)
        
        if isnothing(geo_id)
            error("Unsupported geometry type: $G. Add to geometry_types.")
        end
        write(io, UInt8(geo_id))

        # --- 2. DATA TYPE (T) ---
        type_id = findfirst(==(T), coordinate_types)
        if isnothing(type_id)
            error("Unsupported coordinate type: $T. Add to coordinate_types.")
        end
        write(io, UInt8(type_id))

        # --- 3. REF NODES ---
        ref = m.fe.ref_nodes[1]
        write(io, Int64(length(ref)))
        write(io, ref)

        # --- 4. COORDINATES (x) ---
        write(io, Int64(size(m.x, 1)))
        write(io, Int64(size(m.x, 2)))
        write(io, m.x)

        # --- 5. CONNECTIVITY (el) ---
        # Convert to Int64 for robustness between e.g. 32/64 bit OS
        write(io, Int64(size(m.el, 1)))
        write(io, Int64(size(m.el, 2)))
        write(io, Matrix{Int64}(m.el))

        # --- 6. NEIGHBORS (nb) ---
        write(io, Int64(size(m.nb, 1)))
        write(io, Int64(size(m.nb, 2)))
        write(io, m.nb) 
    end
end

"""
    loadmesh(fname) -> HighOrderMesh

Load a `HighOrderMesh` from a binary `.hom` file previously written by `savemesh`.
"""
function loadmesh(fname)
    open(fname, "r") do io
        # 1. Header Check
        header = String(read(io, 6))
        if header != "HOM_v1"
            error("Invalid file format: $header")
        end

        # 2. Read Metadata
        D = read(io, Int64)
        P = read(io, Int64)
        geo_id = read(io, UInt8)

        # Reconstruct the Geometry Type G (e.g., Simplex{D})
        if geo_id > length(geometry_types)
            error("Unknown geometry ID: $geo_id")
        end
        base_geo = geometry_types[geo_id] # e.g., Simplex
        G = base_geo{D}                   # e.g., Simplex{2}

        # 3. Read Type T
        type_id = read(io, UInt8)
        if type_id > length(coordinate_types)
            error("Unknown coordinate type ID: $type_id")
        end
        T = coordinate_types[type_id]

        # 4. Read ref_nodes
        n_ref = read(io, Int64)
        ref = Vector{T}(undef, n_ref)
        read!(io, ref)

        # 5. Read x
        nx_rows = read(io, Int64)
        nx_cols = read(io, Int64)
        x = Matrix{T}(undef, nx_rows, nx_cols)
        read!(io, x)

        # 6. Read el (Always read as Int64, then cast to system Int)
        nel_rows = read(io, Int64)
        nel_cols = read(io, Int64)
        el_raw = Matrix{Int64}(undef, nel_rows, nel_cols)
        read!(io, el_raw)
        el = Matrix{Int}(el_raw) # Convert to system Int (usually Int64)

        # 7. Read nb
        nnb_rows = read(io, Int64)
        nnb_cols = read(io, Int64)
        nb = Matrix{NeighborData}(undef, nnb_rows, nnb_cols)
        read!(io, nb)

        # Create the HighOrderMesh structure
        fe = FiniteElement(G(), ref)
        HighOrderMesh{D,G,P,T}(fe, x, el, nb)
    end
end

###########################################################################
## ASCII .txt format (legacy, incomplete)
#
# Simple space-delimited text format. Does not store polynomial order,
# reference nodes, or neighbor data — use .hom for full round-trips.

# Write matrix `x` to an open IO stream, one row per line, space-separated.
function write_matrix(f, x)
    for i = 1:size(x,1)
        join(f, x[i,:], ' ')
        println(f)
    end
end

# Read a space-separated matrix from an open IO stream into pre-allocated `x`.
function read_matrix!(f, x)
    for i = 1:size(x,1)
        x[i,:] = parse.(eltype(x), split(readline(f), ' ')[1:size(x,2)])
    end
end

"""
    savemeshtxt(fname, m::HighOrderMesh)

Save the node coordinates and element connectivity of `m` to a plain-text
`.txt` file. Only stores `x` and `el`; polynomial order, reference nodes,
and neighbor data are not saved. Use `savemesh` / `loadmesh` for full round-trips.
"""
function savemeshtxt(fname, m::HighOrderMesh)
    open(fname, "w") do f
        println(f, "$(size(m.x,1)) $(size(m.el,2)) $(size(m.el,1))  # nnodes nelems nnodes_per_elem")
        write_matrix(f, m.x)
        write_matrix(f, m.el')
    end
end

"""
    loadmeshtxt(fname) -> HighOrderMesh

Load a mesh from a plain-text `.txt` file written by `savemeshtxt`.
Returns a linear (`p=1`) `HighOrderMesh`; neighbor data is recomputed from connectivity.
"""
function loadmeshtxt(fname)
    open(fname, "r") do f
        line = split(readline(f), ' ')
        nx, nel, ne = parse.(Int64, line[1:3])
        x  = zeros(Float64, nx, 2)
        el = zeros(Int64, nel, ne)
        read_matrix!(f, x)
        read_matrix!(f, el)
        HighOrderMesh(x, el')
    end
end

