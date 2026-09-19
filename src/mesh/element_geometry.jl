###########################################################################
## Element geometry types

"""
    ElementGeometry{D}

Abstract base type for `D`-dimensional element shapes. The dimension `D` is a
type parameter so that geometry-specific methods can be selected at compile time
via dispatch on the singleton instances `Simplex{D}()` and `Block{D}()`.
"""
abstract type ElementGeometry{D} end

"""
    Simplex{D} <: ElementGeometry{D}

Simplex element: triangle (D=2) or tetrahedron (D=3). Has `D+1` vertices and
`D+1` faces; face `i` is opposite vertex `i`.

One-dimensional meshes always use `Block{1}`; `Simplex{1}` is not supported,
and the sub-geometries of a simplex in dimensions 0 and 1 are `Block{0}` and
`Block{1}`.
"""
struct Simplex{D} <: ElementGeometry{D} end

"""
    Block{D} <: ElementGeometry{D}

Tensor-product (block) element: point (D=0), line segment (D=1),
quadrilateral (D=2), hexahedron (D=3). Has `2^D` vertices and `2D` faces.
Face `2d-1` is the side `ξ_d = 0` and face `2d` is the side `ξ_d = 1`.
"""
struct Block{D} <: ElementGeometry{D} end

# All concrete geometry types, for iteration/registration.
const geometry_types = [Simplex, Block]

###########################################################################
## Basic properties

dim(::ElementGeometry{D}) where {D} = D::Int

name(::Simplex{D}) where {D} = D > 3 ? "$(D)D simplex" :
    ("point", "line", "triangle", "tetrahedron")[D+1]
name(::Block{D}) where {D} = D > 3 ? "$(D)D block" :
    ("point", "line", "quadrilateral", "hexahedron")[D+1]

function Base.show(io::IO, eg::ElementGeometry)
    print(io, "ElementGeometry: $(dim(eg))D $(name(eg))")
end

# The fallback on ElementGeometry acts as a required-method contract:
# any new subtype must define these or it will error at runtime.
nvertices(::ElementGeometry) = error("Not implemented")
nvertices(::Simplex{D}) where {D} = D + 1
nvertices(::Block{D}) where {D} = 2^D

nfaces(::ElementGeometry) = error("Not implemented")
nfaces(::Simplex{D}) where {D} = D + 1
nfaces(::Block{D}) where {D} = 2 * D

nedges(::ElementGeometry) = error("Not implemented")
nedges(::Simplex{D}) where {D} = binomial(D + 1, 2)
nedges(::Block{D}) where {D} = D * 2^(D - 1)

"""
    nnodes(eg::ElementGeometry, p)

Number of Lagrange nodes of polynomial degree `p` on `eg`: `(p+1)^D` for
blocks and `binomial(p+D, D)` for simplices.
"""
nnodes(::Block{D}, p::Integer) where {D} = (p + 1)^D
nnodes(::Simplex{D}, p::Integer) where {D} = binomial(p + D, D)

###########################################################################
## Reference vertices

"""
    vertices(eg::ElementGeometry)

Reference vertex coordinates as an `nvertices × D` integer matrix. Blocks use
lexicographic tensor order on `[0,1]^D` (first coordinate varying fastest);
simplices list the origin followed by the unit vectors.
"""
function vertices(::Block{D}) where {D}
    v = zeros(Int, 2^D, D)
    for (k, i) in enumerate(Iterators.product(ntuple(_ -> 0:1, D)...)), d in 1:D
        v[k, d] = i[d]
    end
    v
end

vertices(::Simplex{D}) where {D} = vcat(zeros(Int, 1, D), Matrix{Int}(I, D, D))

###########################################################################
## Connectivity maps
#
# facemap: local vertex indices for each face, as a (nfv × nf) matrix.
#   Column j gives the vertex indices (within the element) of face j, listed
#   in the vertex order of the face's sub-geometry and oriented so that the
#   outward normal followed by the face's local axes is right-handed.
#
# edgemap: local vertex indices for each edge, as a (2 × nedges) matrix.
#
# plot_face_order: face traversal order for 2D plotting (closed polygon).

facemap(::ElementGeometry) = error("Not implemented")
facemap(::Simplex{2}) = [[2,3] [3,1] [1,2]]
facemap(::Simplex{3}) = [[2,3,4] [1,4,3] [4,1,2] [3,2,1]]
facemap(::Block{1})   = [[1] [2]]
facemap(::Block{2})   = [[3,1] [2,4] [1,2] [4,3]]
facemap(::Block{3})   = [[3,1,7,5] [2,4,6,8] [1,2,5,6]#=
                       =#[4,3,8,7] [3,4,1,2] [5,6,7,8]]

plot_face_order(::ElementGeometry) = error("Not implemented")
plot_face_order(::Simplex{2}) = [1, 2, 3]
plot_face_order(::Block{2})   = [3, 2, 4, 1]

edgemap(::ElementGeometry) = error("Not implemented")
edgemap(::Simplex{2}) = [[2,3] [3,1] [1,2]]
edgemap(::Simplex{3}) = [[1,2] [1,3] [1,4] [2,3] [2,4] [3,4]]
edgemap(::Block{1})   = reshape([1, 2], 2, 1)
edgemap(::Block{2})   = [[1,2] [2,4] [4,3] [3,1]]
edgemap(::Block{3})   = [[1,2] [2,4] [4,3] [3,1] [1,5] [2,6]#=
                       =#[4,8] [3,7] [5,6] [6,8] [8,7] [7,5]]

###########################################################################
## Sub-geometry and utilities

"""
    subgeom(eg::ElementGeometry, d)

Return the `d`-dimensional sub-geometry (face, edge, vertex type) of `eg`.
Sub-geometries of dimension 0 and 1 are always `Block{0}` and `Block{1}`.
"""
subgeom(::Simplex{D}, d::Integer) where {D} = d >= 2 ? Simplex{d}() : Block{d}()
subgeom(::Block{D},   d::Integer) where {D} = Block{d}()

"""Reference-coordinate centroid of the element (on [0,1]^D)."""
midpoint(::Simplex{D}) where {D} = fill(1/(D+1), D)
midpoint(::Block{D})   where {D} = fill(1/2, D)

"""
    find_elgeom(D, nv)

Infer the element geometry from spatial dimension `D` and vertex count `nv`.
Errors if `nv` does not match any known geometry in dimension `D`. In 1D the
result is always `Block{1}`.
"""
function find_elgeom(D, nv)
    nv == nvertices(Block{D}())   && return Block{D}()
    nv == nvertices(Simplex{D}()) && return Simplex{D}()
    error("Cannot determine element geometry for D=$D, nv=$nv")
end

###########################################################################
## Symmetries of the reference element

"""
    symmetries(eg::ElementGeometry) -> Vector{Tuple{Matrix{Int}, Vector{Int}}}

All affine maps `ξ ↦ Aξ + b` that map the reference element onto itself while
permuting its vertices: the hyperoctahedral group (coordinate permutations and
reflections `ξ_d ↦ 1 - ξ_d`) for blocks, and the permutations of the
barycentric coordinates for simplices. Apply a map to a node matrix with
points as rows as `nodes * A' .+ b'`.

A node set that is invariant under all these maps has the same trace on every
face regardless of the face's orientation, which is what a conforming
high-order mesh needs.
"""
function symmetries(::Block{D}) where {D}
    maps = Tuple{Matrix{Int},Vector{Int}}[]
    for σ in _permutations(D), s in Iterators.product(ntuple(_ -> 0:1, D)...)
        A = zeros(Int, D, D)
        b = zeros(Int, D)
        for k in 1:D
            A[k, σ[k]] = 1 - 2s[k]     # ξ'_k = ξ_σ(k)  or  1 - ξ_σ(k)
            b[k] = s[k]
        end
        push!(maps, (A, b))
    end
    maps
end

function symmetries(::Simplex{D}) where {D}
    # Barycentric coordinates λ = (1 - Σξ, ξ_1, …, ξ_D) = L ξ + e_1.
    # Permute them, λ' = P λ, and map back with ξ' = λ'[2:end] = S λ'.
    L  = vcat(-ones(Int, 1, D), Matrix{Int}(I, D, D))
    e1 = vcat(1, zeros(Int, D))
    S  = hcat(zeros(Int, D, 1), Matrix{Int}(I, D, D))
    maps = Tuple{Matrix{Int},Vector{Int}}[]
    for σ in _permutations(D + 1)
        P = Matrix{Int}(I, D + 1, D + 1)[σ, :]
        push!(maps, (S * P * L, S * P * e1))
    end
    maps
end

# All permutations of 1:n as vectors (n! entries; only used for n <= 4).
function _permutations(n::Integer)
    n == 0 && return [Int[]]
    perms = Vector{Int}[]
    for p in _permutations(n - 1), k in 1:n
        push!(perms, insert!(copy(p), k, n))
    end
    perms
end
