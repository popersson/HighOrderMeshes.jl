###########################################################################
## FiniteElement

"""
    FiniteElement{D,G,T}

Lagrange reference element of dimension `D` on geometry `G` with number type
`T`, used for the geometry map of a `HighOrderMesh`. Stores the polynomial
degree `p`, the reference nodes (`nnodes × D`) and the inverse Vandermonde
matrix `coeff`, so that `shapefcns(ξ) = polybasis(ξ) * coeff`.

Reference nodes must form a conforming node set: contain the vertices, be
invariant under all symmetries of the reference element, and be unisolvent;
see [`check_conforming`](@ref).
"""
struct FiniteElement{D,G<:ElementGeometry{D},T}
    p::Int
    nodes::Matrix{T}
    coeff::Matrix{T}
end

###########################################################################
## Constructors

"""
    FiniteElement(eg::ElementGeometry, nodes::AbstractMatrix; check=true)
    FiniteElement(eg::ElementGeometry, p::Int, T=Float64)
    FiniteElement(eg::Block, s1::AbstractVector; check=true)

Construct a reference element on geometry `eg`:

- from an explicit `nnodes × D` node matrix; the degree is inferred from the
  node count and the node set is validated with [`check_conforming`](@ref)
  unless `check=false`;
- with equispaced nodes of degree `p`;
- for blocks, from the tensor product of a 1D node line `s1` of length `p+1`
  (for example `gauss_lobatto01_nodes(p+1)`), validated unless `check=false`.

Mesh elements must pass the validation. Solver elements on deliberately
non-conforming nodes, such as Gauss-Legendre or Gauss-Radau points, are
built with `check=false`; see [`interpolate`](@ref) for moving data between
the two.
"""
function FiniteElement(eg::ElementGeometry{D}, nodes::AbstractMatrix; check::Bool=true) where {D}
    p = _degree_from_nnodes(eg, size(nodes, 1))
    _finite_element(eg, p, nodes, check)
end

FiniteElement(eg::ElementGeometry, p::Integer, ::Type{T}=Float64) where {T} =
    _finite_element(eg, p, T.(equispaced_nodes(eg, p)), false)

FiniteElement(eg::Block{D}, s1::AbstractVector; check::Bool=true) where {D} =
    _finite_element(eg, length(s1) - 1, tensor_nodes(eg, s1), check)

function _finite_element(eg::ElementGeometry{D}, p::Integer, nodes::AbstractMatrix{T}, check::Bool) where {D,T}
    check && check_conforming(eg, nodes)
    V = polybasis(eg, nodes, p)
    size(V, 1) == size(V, 2) ||
        throw(ArgumentError("$(size(nodes,1)) nodes do not match the $(size(V,2)) basis functions of degree $p"))
    coeff = try
        inv(V)
    catch e
        e isa SingularException && throw(ArgumentError("node set is not unisolvent for degree $p"))
        rethrow()
    end
    FiniteElement{D,typeof(eg),T}(Int(p), Matrix{T}(nodes), Matrix{T}(coeff))
end

# Degree whose node count on eg equals ns; throws if there is none.
function _degree_from_nnodes(eg::ElementGeometry{D}, ns::Integer) where {D}
    D == 0 && ns == 1 && return 0
    for p in 1:ns
        nnodes(eg, p) == ns && return p
        nnodes(eg, p) > ns && break
    end
    throw(ArgumentError("$ns nodes do not match any polynomial degree on a $(name(eg))"))
end

###########################################################################
## Accessors

elgeom(::FiniteElement{D,G,T}) where {D,G,T} = G()
dim(::FiniteElement{D}) where {D} = D
porder(fe::FiniteElement) = fe.p
nnodes(fe::FiniteElement) = size(fe.nodes, 1)
name(fe::FiniteElement) = "p=$(fe.p) " * name(elgeom(fe))

function Base.show(io::IO, fe::FiniteElement)
    print(io, "FiniteElement: $(name(fe)), $(nnodes(fe)) nodes")
end

"""
    ref_nodes(fe::FiniteElement)
    ref_nodes(fe::FiniteElement, d)

Reference nodes of `fe` as an `nnodes × D` matrix, or the nodes on the
canonical `d`-dimensional sub-face `ξ_{d+1} = … = ξ_D = 0` as an `n × d`
matrix (`d = 0` gives a `1 × 0` matrix: the vertex at the origin).
"""
ref_nodes(fe::FiniteElement) = fe.nodes

function ref_nodes(fe::FiniteElement{D}, d::Integer) where {D}
    0 <= d <= D || throw(ArgumentError("sub-dimension $d is outside 0:$D"))
    d == D ? fe.nodes : _restrict_nodes(fe.nodes, d)
end

"""
    subelement(fe::FiniteElement, d)

The `d`-dimensional reference element on `subgeom(elgeom(fe), d)` whose nodes
are the restriction of the nodes of `fe` to the canonical `d`-face.
"""
function subelement(fe::FiniteElement{D}, d::Integer) where {D}
    d == D && return fe
    _finite_element(subgeom(elgeom(fe), d), fe.p, ref_nodes(fe, d), false)
end

"""Indices of the reference vertices within the node set of `fe`."""
function corner_nodes(fe::FiniteElement)
    v  = vertices(elgeom(fe))
    ix = [ _find_node(fe.nodes, view(v, i, :)) for i in axes(v, 1) ]
    any(isnothing, ix) && error("reference vertices are missing from the node set")
    Int.(ix)
end

# Tolerance for coordinate comparisons: exact types compare exactly.
_tol(::Type{T}) where {T<:AbstractFloat} = sqrt(eps(T))
_tol(::Type) = 0

# Index of the node coinciding with the point x (rows of nodes), or nothing.
function _find_node(nodes::AbstractMatrix{T}, x) where {T}
    tol = _tol(T)
    for k in axes(nodes, 1)
        all(abs(nodes[k, d] - x[d]) <= tol for d in axes(nodes, 2)) && return k
    end
    nothing
end

# Rows of nodes lying on the canonical d-face, keeping the first d coordinates.
function _restrict_nodes(nodes::AbstractMatrix{T}, d::Integer) where {T}
    tol  = _tol(T)
    keep = [ all(abs(nodes[k, e]) <= tol for e in d+1:size(nodes, 2)) for k in axes(nodes, 1) ]
    nodes[keep, 1:d]
end

###########################################################################
## Conformity check

"""
    check_conforming(eg::ElementGeometry, nodes::AbstractMatrix)

Verify that `nodes` (`nnodes × D`) can serve as the reference node set of a
conforming high-order mesh on `eg`, and throw an `ArgumentError` otherwise.
The node set must

1. have the node count of some polynomial degree `p ≥ 1`,
2. contain the reference vertices,
3. be invariant under every map in `symmetries(eg)`, so that its trace on a
   face does not depend on the orientation of the face,
4. be unisolvent for degree `p`, and
5. restrict to a conforming node set of the sub-geometry on the canonical
   face `ξ_D = 0` (checked recursively down to the edges).
"""
function check_conforming(eg::ElementGeometry{D}, nodes::AbstractMatrix{T}) where {D,T}
    size(nodes, 2) == D ||
        throw(ArgumentError("node matrix has $(size(nodes, 2)) columns, expected D=$D"))
    p = _degree_from_nnodes(eg, size(nodes, 1))

    for (i, v) in enumerate(eachrow(vertices(eg)))
        isnothing(_find_node(nodes, v)) &&
            throw(ArgumentError("node set does not contain reference vertex $i at $(collect(v))"))
    end

    for (A, b) in symmetries(eg)
        mapped = nodes * A' .+ b'
        for k in axes(mapped, 1)
            isnothing(_find_node(nodes, view(mapped, k, :))) &&
                throw(ArgumentError("node set is not symmetric: the image of node $k under a symmetry of the element is not a node"))
        end
    end

    V = polybasis(eg, nodes, p)
    issuccess(lu(V; check=false)) ||
        throw(ArgumentError("node set is not unisolvent for degree $p"))
    T <: AbstractFloat && cond(V) > 1e12 &&
        throw(ArgumentError("node set is numerically not unisolvent for degree $p (condition number $(cond(V)))"))

    if D >= 1
        sg  = subgeom(eg, D - 1)
        sub = _restrict_nodes(nodes, D - 1)
        size(sub, 1) == nnodes(sg, p) ||
            throw(ArgumentError("the face ξ_$D = 0 carries $(size(sub, 1)) nodes; a degree-$p $(name(sg)) needs $(nnodes(sg, p))"))
        check_conforming(sg, sub)
    end
    nothing
end

###########################################################################
## Shape functions and interpolation

"""
    shapefcns(fe::FiniteElement, ξ)
    dshapefcns(fe::FiniteElement, ξ)
    shapefcns(eg::ElementGeometry, ξ)
    dshapefcns(eg::ElementGeometry, ξ)

Lagrange shape functions of `fe` at the reference points `ξ` (`nξ × D`), as
an `nξ × nnodes` matrix, and their reference gradients as an
`nξ × nnodes × D` array. The forms taking a geometry use its linear (`p=1`)
element.
"""
shapefcns(fe::FiniteElement, ξ::AbstractVecOrMat) = polybasis(elgeom(fe), ξ, fe.p) * fe.coeff

function dshapefcns(fe::FiniteElement{D}, ξ::AbstractVecOrMat) where {D}
    dV  = dpolybasis(elgeom(fe), ξ, fe.p)
    out = similar(dV, promote_type(eltype(dV), eltype(fe.coeff)), (size(dV, 1), nnodes(fe), D))
    for d in 1:D
        out[:, :, d] = dV[:, :, d] * fe.coeff
    end
    out
end

shapefcns(eg::ElementGeometry, ξ::AbstractVecOrMat)  = shapefcns(_linear_element(eg, ξ), ξ)
dshapefcns(eg::ElementGeometry, ξ::AbstractVecOrMat) = dshapefcns(_linear_element(eg, ξ), ξ)

_linear_element(eg, ξ) = FiniteElement(eg, 1, eltype(ξ) <: Integer ? Float64 : eltype(ξ))

"""
    interpolate(N, u)
    interpolate(dN, u)

Apply the shape-function matrix `N = shapefcns(fe, ξ)` (`nξ × ns`) to nodal
values `u` (`ns × …`, for example `ns × nel` or `ns × nel × ncomp`), giving
the interpolated values at the points `ξ` as `nξ × …`. With the gradient
array `dN = dshapefcns(fe, ξ)` (`nξ × ns × D`) the result has an extra
trailing dimension `D` holding the reference derivatives.

```julia
N  = shapefcns(m.fe, ξ)
x  = interpolate(N, dg_nodes(m))                 # nξ × nel × D physical coordinates
J  = interpolate(dshapefcns(m.fe, ξ), dg_nodes(m))   # nξ × nel × D × D, J[:,:,i,j] = ∂x_i/∂ξ_j
```
"""
function interpolate(N::AbstractMatrix, u::AbstractArray)
    ns = size(u, 1)
    size(N, 2) == ns || throw(DimensionMismatch("N has $(size(N, 2)) columns but u has $ns rows"))
    reshape(N * reshape(u, ns, :), size(N, 1), Base.tail(size(u))...)
end

function interpolate(dN::AbstractArray{<:Any,3}, u::AbstractArray)
    D   = size(dN, 3)
    out = similar(u, promote_type(eltype(dN), eltype(u)), (size(dN, 1), Base.tail(size(u))..., D))
    for d in 1:D
        selectdim(out, ndims(out), d) .= interpolate(view(dN, :, :, d), u)
    end
    out
end

"""
    interpolate(fe_from::FiniteElement, fe_to::FiniteElement, u)

Re-nodalize nodal values `u` (`nnodes(fe_from) × …`) from the node set of
`fe_from` to the node set of `fe_to`, both on the same reference geometry.
When the two elements have the same degree the result is exact, because the
Lagrange interpolants through either node set are the same polynomial; for a
lower-degree `fe_to` it is the interpolation of the field at the coarser
nodes.

Typical use: a solver working on its own non-conforming nodes transfers its
solution to the mesh nodes for plotting or export.

```julia
fe_sol = FiniteElement(Block{2}(), gauss_legendre01_nodes(p+1); check=false)
u_mesh = interpolate(fe_sol, m.fe, u_sol)    # nnodes(m.fe) × nel
```
"""
function interpolate(fe_from::FiniteElement{D}, fe_to::FiniteElement{D}, u::AbstractArray) where {D}
    elgeom(fe_from) == elgeom(fe_to) ||
        throw(ArgumentError("cannot interpolate between a $(name(elgeom(fe_from))) and a $(name(elgeom(fe_to)))"))
    interpolate(shapefcns(fe_from, ref_nodes(fe_to)), u)
end
