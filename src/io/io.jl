#######################################################
## File I/O - binary .hom format

using Serialization

"""
    savemesh(fname, m::HighOrderMesh{D,G,P,T})

Save HighOrderMesh to .hom format.
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
        serialize(io, T)

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
    loadmesh(fname) -> HighOrderMesh{D,G,P,T}

Load HighOrderMesh from .hom format.
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
        T = deserialize(io)

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

#######################################################
## File I/O - simple but incomplete .txt format

"""
    write_matrix(f, x)

TBW
"""
function write_matrix(f, x)
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            print(f, x[i,j])
            if j<size(x,2)
                print(f, ' ')
            end
        end
        println(f)
    end
end

"""
    read_matrix!(f, x)

TBW
"""
function read_matrix!(f, x)
    for i = 1:size(x,1)
        line = readline(f)
        x[i,:] = parse.(eltype(x), split(line, ' ')[1:size(x,2)])
    end
end

"""
    savemeshtxt(fname, m::HighOrderMesh)

TBW
"""
function savemeshtxt(fname, m::HighOrderMesh)
    nx = size(m.x,1)
    nel = size(m.el,2)

    f = open(fname, "w")
    println(f, "$nx $nel $(size(m.el,1))  # nbr_nodes nbr_elems nbr_nodes_per_elem")
    write_matrix(f, m.x)
    write_matrix(f, m.el')
    close(f)
end

"""
    loadmeshtxt(fname)

TBW
"""
function loadmeshtxt(fname)
    f = open(fname, "r")
    line = split(readline(f), " ")
    nx,nel,ne = parse.(Int64, line[1:3])
    name = line[3]
    x = zeros(Float64, nx, 2)
    read_matrix!(f, x)
    el = zeros(Int64, nel, ne)
    read_matrix!(f, el)
    close(f)
    return HighOrderMesh(x,el')
end

