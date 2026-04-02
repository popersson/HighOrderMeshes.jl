# (neighbor element, neighbor face index, orientation flag)
const NeighborData = Tuple{Int32, Int16, Int16}

"""
    HighOrderMesh{D,G,P,T}

High-order unstructured mesh of dimension `D`, element geometry `G`,
polynomial order `P`, and floating-point type `T`.

- `fe`: reference element (basis and quadrature nodes)
- `x`:  global node coordinates (`nnodes × D`)
- `el`: element-to-node connectivity (`nnodes_per_elem × nelems`)
- `nb`: neighbor data per face (`nfaces × nelems`); boundary faces have
        a non-positive first entry (`-bnd_number`, 0, 0).
"""
struct HighOrderMesh{D,G,P,T}
    fe::FiniteElement{D,G,P,T}
    x::Matrix{T}
    el::Matrix{Int}
    nb::Matrix{NeighborData}
end

###########################################################################
## Constructors

"""
    HighOrderMesh(fe, x, el; bndexpr)
    HighOrderMesh(x, el; bndexpr)

Construct a `HighOrderMesh` from node coordinates `x` and element connectivity `el`.

- First form: supply a `FiniteElement` `fe` for an arbitrary polynomial order.
- Second form: geometry and order are inferred from `x` and `el` (linear, `p=1`).

`bndexpr` is a function `x -> [expr1(x), expr2(x), ...]` where each `expri`
evaluates to zero on boundary region `i`. Defaults to a single region.
"""
function HighOrderMesh(fe::FiniteElement{D,G,P,T},
                       x::AbstractMatrix{T},
                       el::AbstractMatrix{Int};
                       bndexpr=p->[0]) where {D,G,P,T}
    nb = el2nb(el[corner_nodes(fe),:], G())
    m = HighOrderMesh{D,G,P,T}(fe, x, el, nb)
    set_bnd_numbers!(m, bndexpr)
    m
end

function HighOrderMesh(x::Matrix{T}, el::AbstractMatrix{Int}; kwargs...) where {T}
    dim, nv = size(x,2), size(el,1)
    eg = find_elgeom(dim, nv)
    fe = FiniteElement(eg, 1, T)
    HighOrderMesh(fe, x, el; kwargs...)
end

###########################################################################
## Accessors

"""Element-local DG node coordinates: `nnodes_per_elem × nelems × D` array."""
dg_nodes(m::HighOrderMesh) = m.x[m.el,:]

elgeom(::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = G()
dim(::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = D
porder(::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = P
nnodes(m::HighOrderMesh) = size(m.x,1)
nel(m::HighOrderMesh) = size(m.el,2)

function Base.show(io::IO, m::HighOrderMesh)
    print(io, "HighOrderMesh: $(dim(m))D, ")
    print(io, "$(nnodes(m)) vertices, ")
    print(io, "$(nel(m)) $(name(m.fe)) elements.")
end

###########################################################################
## Neighbor connectivity

# Build the neighbor matrix from the linear corner-node connectivity.
# Faces are matched by sorting their vertex indices; unmatched faces are
# boundary faces and remain (0,0,0) until set_bnd_numbers! labels them.
# TODO: implement neighbor face permutation (orientation tracking).
function el2nb(el, eg)
    fmap = facemap(eg)
    nv, nel = size(el)
    nfv, nf = size(fmap)

    nb = fill(NeighborData((0,0,0)), nf, nel)
    dd = Dict{NTuple{nfv,Int}, NeighborData}()
    sizehint!(dd, nel * nf)
    e = fill(0, nfv)
    for iel = 1:nel
        for jf = 1:nf
            e[:] = el[fmap[:,jf], iel]
            sort!(e)
            et = Tuple(e)
            if haskey(dd, et)
                nbel = pop!(dd, et)
                nb[jf,iel]          = nbel
                nb[nbel[2],nbel[1]] = (iel, jf, 0)
            else
                dd[et] = (iel, jf, 0)
            end
        end
    end
    nb
end

###########################################################################
## Polynomial order and node set changes

"""
    set_ref_nodes(m::HighOrderMesh, newfe::FiniteElement)
    set_ref_nodes(m::HighOrderMesh{D,G,P,T}, newsline::Vector{T})
    set_degree(m::HighOrderMesh, newp::Int)
    set_lobatto_nodes(m::HighOrderMesh{D,Block{D},P,T})

Return a new mesh with the reference nodes (and polynomial order) changed.

`set_ref_nodes` is the general form; `set_degree` and `set_lobatto_nodes` are
convenience wrappers. Physical node coordinates are recomputed by evaluating
the old basis at the new reference nodes and deduplicating.
"""
function set_ref_nodes(m::HighOrderMesh{D,G,P,T}, newfe::FiniteElement) where {D,G,P,T}
    newp  = porder(newfe)
    newns = nbr_ho_nodes(newfe)
    nv, nel = size(m.el)
    # Evaluate old basis at new reference nodes to get new physical coordinates
    Pfe    = eval_poly(G(), ref_nodes(newfe,D), P)
    xdg    = dg_nodes(m)
    newxdg = reshape(Pfe * (m.fe.shapefcn_coeff[D] * reshape(xdg, nv, D*nel)), newns*nel, D)
    neweldg = reshape(1:newns*nel, newns, nel)
    newx, newel = unique_mesh_nodes(newxdg, neweldg)
    HighOrderMesh{D,G,newp,T}(newfe, newx, newel, m.nb)
end

set_ref_nodes(m::HighOrderMesh{D,G,P,T}, newsline::Vector{T}) where {D,G,P,T} =
    set_ref_nodes(m, FiniteElement(G(), newsline))

"""Change polynomial order, keeping equispaced nodes."""
set_degree(m::HighOrderMesh{D,G,P,T}, newp::Int) where {D,G,P,T} =
    set_ref_nodes(m, FiniteElement(G(), newp, T))

"""Switch to Gauss-Lobatto nodes (Block elements only)."""
set_lobatto_nodes(m::HighOrderMesh{D,Block{D},P,T}) where {D,P,T} =
    set_ref_nodes(m, gauss_lobatto01_nodes(P+1))

###########################################################################
## Face-to-node mapping

"""
    mkface2nodes(eg::ElementGeometry, sface, svol)
    mkface2nodes(fe::FiniteElement)
    mkface2nodes(m::HighOrderMesh)

Build the face-to-node index map for element geometry `eg`. Returns an
`(nfacenodes × nfaces)` integer matrix where column `j` contains the local
node indices (within the element) that lie on face `j`.

Matching is done by comparing shape function values at face reference nodes
against volume reference nodes, so it is basis-order-aware.
"""
function mkface2nodes(eg::ElementGeometry{D}, sface::AbstractArray{T}, svol::AbstractArray{T}) where {D,T}
    fmap       = facemap(eg)
    basis_face = snap.(eval_shapefcns(subgeom(eg, D-1), sface))
    basis_vol  = snap.(eval_shapefcns(eg, svol))
    f2n = fill(0, size(basis_face,1), size(fmap,2))
    for (ii, ic) in enumerate(eachcol(fmap))
        f2n[:,ii] .= indexin(eachrow(basis_face), eachrow(basis_vol[:,ic]))
    end
    f2n
end

mkface2nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T} =
    mkface2nodes(G(), ref_nodes(fe, D-1), ref_nodes(fe, D))

mkface2nodes(m::HighOrderMesh) = mkface2nodes(m.fe)
    