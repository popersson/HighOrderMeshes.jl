# (neighbor element, neighbor face index, one-based face permutation)
const NeighborData = Tuple{Int32, Int16, Int16}

"""
    HighOrderMesh{D,G,T}

High-order unstructured mesh of dimension `D`, element geometry `G`, and
number type `T`. All elements share one reference element `fe`, so the
mesh has a single element type and a single polynomial degree.

- `fe`: reference element (degree, reference nodes, shape functions)
- `x`:  global node coordinates (`nnodes × D`)
- `el`: element-to-node connectivity (`nnodes_per_elem × nelems`)
- `nb`: neighbor data per face (`nfaces × nelems`). Interior faces store
        `(neighbor element, neighbor face, one-based face permutation)`;
        boundary faces store `(-bnd_number, 0, 0)`.

The element-local (DG) node coordinates are `x[el, :]`, see [`dg_nodes`](@ref).
"""
struct HighOrderMesh{D,G<:ElementGeometry{D},T}
    fe::FiniteElement{D,G,T}
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

- First form: supply a `FiniteElement` `fe` for an arbitrary polynomial degree.
- Second form: geometry and degree are inferred from `x` and `el` (linear, `p=1`).

`bndexpr` is a function `x -> [expr1(x), expr2(x), ...]` where each `expri`
evaluates to zero on boundary region `i`. Defaults to a single region.
"""
function HighOrderMesh(fe::FiniteElement{D,G,T},
                       x::AbstractMatrix{T},
                       el::AbstractMatrix{Int};
                       bndexpr=p->[0]) where {D,G,T}
    nb = el2nb(el[corner_nodes(fe),:], G())
    m = HighOrderMesh{D,G,T}(fe, x, el, nb)
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

elgeom(::HighOrderMesh{D,G,T}) where {D,G,T} = G()
dim(::HighOrderMesh{D}) where {D} = D
porder(m::HighOrderMesh) = porder(m.fe)
nnodes(m::HighOrderMesh) = size(m.x,1)
nel(m::HighOrderMesh) = size(m.el,2)

function Base.show(io::IO, m::HighOrderMesh)
    print(io, "HighOrderMesh: $(dim(m))D, ")
    print(io, "$(nnodes(m)) nodes, ")
    print(io, "$(nel(m)) $(name(m.fe)) elements.")
end

###########################################################################
## Neighbor connectivity

# One-based face permutation matching 3DG's convention: for a shared 3D face,
# store the position of this face's first node in the neighbor's face ordering.
function face_permutation(::ElementGeometry{D}, f1, f2) where {D}
    D <= 2 && return Int16(1)
    perm = findfirst(==(first(f1)), f2)
    isnothing(perm) && error("Could not determine neighbor face permutation")
    Int16(perm)
end

# Build the neighbor matrix from the linear corner-node connectivity.
# Faces are matched by sorting their vertex indices; unmatched faces are
# boundary faces and remain (0,0,0) until set_bnd_numbers! labels them.
function el2nb(el, eg)
    fmap = facemap(eg)
    nv, nel = size(el)
    nfv, nf = size(fmap)

    nb = fill(NeighborData((0,0,0)), nf, nel)
    dd = Dict{NTuple{nfv,Int}, NeighborData}()
    sizehint!(dd, nel * nf)
    e = fill(0, nfv)
    f = fill(0, nfv)
    for iel = 1:nel
        for jf = 1:nf
            f[:] = el[fmap[:,jf], iel]
            e[:] = f
            et = Tuple(sort!(e))
            if haskey(dd, et)
                nbel = pop!(dd, et)
                f2 = el[fmap[:,nbel[2]], nbel[1]]
                nb[jf,iel]          = (nbel[1], nbel[2], face_permutation(eg, f, f2))
                nb[nbel[2],nbel[1]] = (iel, jf, face_permutation(eg, f2, f))
            else
                dd[et] = (iel, jf, 0)
            end
        end
    end
    nb
end

###########################################################################
## Polynomial degree and node set changes

"""
    set_ref_nodes(m::HighOrderMesh, newfe::FiniteElement)
    set_ref_nodes(m::HighOrderMesh, nodes::AbstractMatrix)
    set_ref_nodes(m::HighOrderMesh{D,Block{D}}, s1::AbstractVector)
    set_degree(m::HighOrderMesh, newp::Int)
    set_lobatto_nodes(m::HighOrderMesh{D,Block{D}})

Return a new mesh with the reference nodes (and polynomial degree) changed.
The new mesh shares no arrays with `m`.

`set_ref_nodes` is the general form, taking a reference element, an explicit
`nnodes × D` node matrix, or (for blocks) a symmetric 1D node line.
`set_degree` switches to equispaced nodes of degree `newp` and
`set_lobatto_nodes` to the tensor product of Gauss-Lobatto nodes. Physical
node coordinates are recomputed by evaluating the old geometry map at the new
reference nodes and deduplicating.
"""
function set_ref_nodes(m::HighOrderMesh{D,G,T}, newfe::FiniteElement{D,G,T}) where {D,G,T}
    N       = shapefcns(m.fe, ref_nodes(newfe))       # nnew × nold
    newxdg  = interpolate(N, dg_nodes(m))             # nnew × nel × D
    newns   = nnodes(newfe)
    neweldg = reshape(1:newns*nel(m), newns, nel(m))
    newx, newel = unique_mesh_nodes(reshape(newxdg, newns*nel(m), D), neweldg)
    HighOrderMesh{D,G,T}(newfe, newx, newel, copy(m.nb))
end

set_ref_nodes(m::HighOrderMesh{D,G,T}, nodes::AbstractMatrix) where {D,G,T} =
    set_ref_nodes(m, FiniteElement(G(), T.(nodes)))

set_ref_nodes(m::HighOrderMesh{D,Block{D},T}, s1::AbstractVector) where {D,T} =
    set_ref_nodes(m, FiniteElement(Block{D}(), T.(s1)))

"""Change polynomial degree, keeping equispaced nodes."""
set_degree(m::HighOrderMesh{D,G,T}, newp::Integer) where {D,G,T} =
    set_ref_nodes(m, FiniteElement(G(), newp, T))

"""Switch to Gauss-Lobatto nodes (Block elements only)."""
set_lobatto_nodes(m::HighOrderMesh{D,Block{D},T}) where {D,T} =
    set_ref_nodes(m, gauss_lobatto01_nodes(porder(m)+1, T))

###########################################################################
## Face-to-node mapping

"""
    mkface2nodes(eg::ElementGeometry, sface, svol)
    mkface2nodes(fe::FiniteElement)
    mkface2nodes(m::HighOrderMesh)

Build the face-to-node index map for element geometry `eg`. Returns an
`(nfacenodes × nfaces)` integer matrix where column `j` contains the local
node indices (within the element) that lie on face `j`.

Matching is done by comparing linear shape function values at the face
reference nodes `sface` against the volume reference nodes `svol`, so it is
independent of the node order.
"""
function mkface2nodes(eg::ElementGeometry{D}, sface::AbstractMatrix, svol::AbstractMatrix) where {D}
    fmap       = facemap(eg)
    basis_face = snap.(shapefcns(subgeom(eg, D-1), sface))
    basis_vol  = snap.(shapefcns(eg, svol))
    f2n = fill(0, size(basis_face,1), size(fmap,2))
    for (ii, ic) in enumerate(eachcol(fmap))
        f2n[:,ii] .= indexin(eachrow(basis_face), eachrow(basis_vol[:,ic]))
    end
    f2n
end

mkface2nodes(fe::FiniteElement{D}) where {D} =
    mkface2nodes(elgeom(fe), ref_nodes(fe, D-1), ref_nodes(fe))

mkface2nodes(m::HighOrderMesh) = mkface2nodes(m.fe)
