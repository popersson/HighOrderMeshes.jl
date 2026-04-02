###########################################################################
## FiniteElement

"""
    FiniteElement{D,G,P,T}

Reference element of dimension `D`, geometry `G`, polynomial order `P`, and
floating-point type `T`. Stores tensor-product reference nodes and the
pre-factored inverse Vandermonde matrices used for shape function evaluation.
"""
struct FiniteElement{D,G<:ElementGeometry{D},P,T}
    ref_nodes::NTuple{D, Matrix{T}}      # ref_nodes[d]: nodes on d-dim sub-geometry
    shapefcn_coeff::NTuple{D, Matrix{T}} # inv(Vandermonde), one per dimension
end

# Accessors
elgeom(::FiniteElement{D,G,P,T}) where {D,G,P,T} = G()
dim(::FiniteElement{D,G,P,T}) where {D,G,P,T} = D
porder(::FiniteElement{D,G,P,T}) where {D,G,P,T} = P
name(::FiniteElement{D,G,P,T}) where {D,G,P,T} = "p=$(P) " * name(G())

"""Return reference nodes for the `dim`-dimensional sub-geometry."""
ref_nodes(fe::FiniteElement{D,G,P,T}, dim) where {D,G,P,T} =
    dim == 0 ? zeros(T,1,0) : fe.ref_nodes[dim]

"""Number of nodes on the `dims`-dimensional sub-geometry (defaults to full element)."""
nbr_ho_nodes(fe::FiniteElement{D}, dims=D) where {D} = size(ref_nodes(fe,dims),1)

"""Indices of the linear corner nodes within the high-order node set."""
function corner_nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T}
    fe1 = FiniteElement(G(), 1)
    indexin(eachrow(ref_nodes(fe1,D)), eachrow(ref_nodes(fe,D)))
end

function Base.show(io::IO, fe::FiniteElement)
    print(io, "FiniteElement: $(dim(fe))D, $(name(fe)) element")
end

###########################################################################
## Constructors

"""
    FiniteElement(eg::ElementGeometry{D}, sline::AbstractVector{T})
    FiniteElement(eg::ElementGeometry, p::Int, T=Float64)

Construct a `FiniteElement` on geometry `eg`.

- First form: use the 1D node distribution `sline` (length `p+1`) to build
  tensor-product reference nodes and invert the Vandermonde matrix.
- Second form: use equispaced nodes of polynomial order `p`.
"""
function FiniteElement(eg::ElementGeometry{D}, sline::AbstractVector{T}) where {D,T}
    p = length(sline) - 1
    nodes  = Tuple( ref_nodes(subgeom(eg,d), sline) for d in 1:D )
    coeffs = Tuple( inv(eval_poly(eg, nodes[d], p)) for d in 1:D )
    FiniteElement{D,typeof(eg),p,T}(nodes, coeffs)
end

FiniteElement(eg::ElementGeometry, p::Int, T=Float64) = FiniteElement(eg, T.(equispaced(p)))

# Evaluate the polynomial basis for this element at points s.
eval_poly(::FiniteElement{D,G,P,T}, s; gradient=false) where {D,G,P,T} =
    eval_poly(G(), s, P; gradient=gradient)

###########################################################################
## Shape functions

"""
    eval_shapefcns(fe::FiniteElement, ss; gradient=false)
    eval_shapefcns(eg::ElementGeometry, ss; gradient=false)

Evaluate shape functions at reference coordinates `ss` (`nss × ndim`).

- `gradient=false` → `nss × ns`
- `gradient=true`  → `nss × ns × ndim`

The second form uses a linear (`p=1`) element on geometry `eg`.
"""
function eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    nss, ndim = size(ss,1), size(ss,2)
    ndim == 0 && return gradient ? ones(T,nss,1,0) : ones(T,nss,1)
    pol = eval_poly(G(), ss, P; gradient=gradient)
    C = fe.shapefcn_coeff[ndim]
    return gradient ? cat( (pol[:,:,k] * C for k = 1:ndim)..., dims=3) : pol * C
end

eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}; gradient=false) where {T} =
    eval_shapefcns(FiniteElement(eg, 1, T), ss; gradient=gradient)

"""
    eval_field(fe::FiniteElement, u, ss; gradient=false)

Evaluate the FEM field `u` at reference coordinates `ss` for all elements.

- `u`:  `ns × nel` or `ns × nel × nc` (nodal DOFs, optionally `nc` components)
- `ss`: `nss × ndim` (reference coordinates)
- `gradient=false` → `nss × nel` or `nss × nel × nc`
- `gradient=true`  → `nss × nel × nc × ndim`
"""
function eval_field(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    nss, ndim = size(ss,1), size(ss,2)
    ns, nel   = size(u,1), size(u,2)
    @assert ns == nbr_ho_nodes(fe, ndim)

    ndim == 0 && return gradient ? zeros(T,nss,nel,size(u,3),0) : repeat(u[[1],:,:],nss)

    C = fe.shapefcn_coeff[ndim]
    V = eval_poly(fe, ss; gradient=gradient)

    if gradient
        V = reshape(permutedims(V,(1,3,2)), :, ns)  # (nss*ndim) × ns
    end
    matrix_output = V * (C * reshape(u, ns, :))     # apply basis and DOFs together

    if gradient
        return permutedims(reshape(matrix_output, nss, ndim, nel, size(u,3)), (1,3,4,2))
    elseif ndims(u) <= 2
        return matrix_output
    else
        return reshape(matrix_output, nss, nel, :)
    end
end

#################################################################################
