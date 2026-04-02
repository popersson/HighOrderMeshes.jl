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

Simplex element: point (D=0), line segment (D=1), triangle (D=2),
tetrahedron (D=3). Has `D+1` vertices and `D+1` faces.
"""
struct Simplex{D} <: ElementGeometry{D} end

"""
    Block{D} <: ElementGeometry{D}

Tensor-product (block) element: point (D=0), line segment (D=1),
quadrilateral (D=2), hexahedron (D=3). Has `2^D` vertices and `2D` faces.
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

###########################################################################
## Connectivity maps
#
# facemap: local vertex indices for each face, as a (nfv × nf) matrix.
#   Column j gives the vertex indices (within the element) of face j.
#   Ordering follows the outward-normal convention for each geometry.
#
# edgemap: local vertex indices for each edge, as a (2 × nedges) matrix.
#
# plot_face_order: face traversal order for 2D plotting (closed polygon).

facemap(::ElementGeometry) = error("Not implemented")
facemap(::Simplex{1}) = [[2] [1]]
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
edgemap(::Simplex{1}) = [[2,1]]
edgemap(::Simplex{2}) = [[2,3] [3,1] [1,2]]
edgemap(::Simplex{3}) = [[1,2] [1,3] [1,4] [2,3] [2,4] [3,4]]
edgemap(::Block{1})   = [[2,1]]
edgemap(::Block{2})   = [[1,2] [2,4] [4,3] [3,1]]
edgemap(::Block{3})   = [[1,2] [2,4] [4,3] [3,1] [1,5] [2,6]#=
                       =#[4,8] [3,7] [5,6] [6,8] [8,7] [7,5]]

###########################################################################
## Sub-geometry and utilities

"""Return the `newD`-dimensional sub-geometry (face/edge type) of `eg`."""
subgeom(::Simplex{D}, newD) where {D} = Simplex{newD}()
subgeom(::Block{D},   newD) where {D} = Block{newD}()

"""Reference-coordinate centroid of the element (on [0,1]^D)."""
midpoint(::Simplex{D}) where {D} = fill(1/(D+1), D)
midpoint(::Block{D})   where {D} = fill(1/2, D)

"""
    find_elgeom(D, nv)

Infer the element geometry from spatial dimension `D` and vertex count `nv`.
Errors if `nv` does not match any known geometry in dimension `D`.
"""
function find_elgeom(D, nv)
    nv == nvertices(Block{D}())   && return Block{D}()
    nv == nvertices(Simplex{D}()) && return Simplex{D}()
    error("Cannot determine element geometry for D=$D, nv=$nv")
end
