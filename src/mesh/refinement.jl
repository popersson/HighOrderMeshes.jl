###########################################################################
## Isoparametric h-refinement of 2D meshes
#
# An element is refined by a pattern chosen from its marked faces. A pattern
# lists the child elements by their vertices in the parent's reference
# coordinates, and each child inherits the parent's geometry by evaluating
# the parent's shape functions at the child's reference nodes. The quad
# patterns and the child order follow 3DG's qmshrefine: the first child takes
# the parent's element index and the others are appended to the mesh, in
# element order.

# Marked faces => (pattern vertices in units of 1/2, children as rows of
# vertex indices in the vertex order of the element geometry).
_refinement_patterns(::Block{2}) = [
    Bool[0,0,1,1] => ([0 0; 2 0; 0 2; 2 2; 1 0; 1 2], [5 2 6 4; 1 5 3 6]),
    Bool[1,1,0,0] => ([0 0; 2 0; 0 2; 2 2; 0 1; 2 1], [5 6 3 4; 1 2 5 6]),
    Bool[1,0,1,0] => ([0 0; 2 0; 0 2; 2 2; 1 0; 0 1; 1 1], [6 7 3 4; 1 5 6 7; 5 2 7 4]),
    Bool[1,0,0,1] => ([0 0; 2 0; 0 2; 2 2; 0 1; 1 2; 1 1], [7 2 6 4; 1 2 5 7; 5 7 3 6]),
    Bool[0,1,1,0] => ([0 0; 2 0; 0 2; 2 2; 1 0; 2 1; 1 1], [7 6 3 4; 1 5 3 7; 5 2 7 6]),
    Bool[0,1,0,1] => ([0 0; 2 0; 0 2; 2 2; 2 1; 1 2; 1 1], [7 5 6 4; 1 2 7 5; 1 7 3 6]),
    Bool[1,1,1,1] => ([0 0; 2 0; 0 2; 2 2; 1 0; 0 1; 2 1; 1 2; 1 1],
                      [1 5 6 9; 5 2 9 7; 6 9 3 8; 9 7 8 4]),
]

_refinement_patterns(::Simplex{2}) = [
    Bool[1,0,0] => ([0 0; 2 0; 0 2; 1 1], [1 2 4; 1 4 3]),
    Bool[0,1,0] => ([0 0; 2 0; 0 2; 0 1], [2 3 4; 2 4 1]),
    Bool[0,0,1] => ([0 0; 2 0; 0 2; 1 0], [3 1 4; 3 4 2]),
    Bool[1,1,1] => ([0 0; 2 0; 0 2; 1 1; 0 1; 1 0], [1 6 5; 6 2 4; 5 4 3; 4 5 6]),
]

# The parent face containing the points ξ (rows), or 0 if there is none:
# the linear shape functions of the vertices off the face vanish on it.
function _parent_face(eg::ElementGeometry, ξ)
    fmap = facemap(eg)
    N    = shapefcns(eg, ξ)
    tol  = sqrt(eps(float(eltype(N))))
    onface(j) = all(abs.(N[:, setdiff(1:nvertices(eg), fmap[:,j])]) .< tol)
    something(findfirst(onface, axes(fmap,2)), 0)
end

# For each pattern: the children's interpolation matrices from the parent
# nodes, and for each child face the parent face it lies on (0 if interior).
function _refinement_rules(fe::FiniteElement{2,G,T}) where {G,T}
    eg    = G()
    fmap  = facemap(eg)
    ξlin  = shapefcns(eg, ref_nodes(fe))          # ns × nv
    rules = Dict{Vector{Bool}, @NamedTuple{N::Vector{Matrix{T}}, faces::Matrix{Int}}}()
    for (marks, (pts, children)) in _refinement_patterns(eg)
        v     = T.(pts) ./ 2
        N     = [ shapefcns(fe, ξlin * v[c,:]) for c in eachrow(children) ]
        faces = [ _parent_face(eg, v[c[fmap[:,j]],:]) for j in axes(fmap,2), c in eachrow(children) ]
        rules[marks] = (; N, faces)
    end
    rules
end

# Mark both sides of every interior face that is marked on either side.
function _share_marks!(marked, m::HighOrderMesh)
    for iel in axes(marked,2), j in axes(marked,1)
        nb = m.nb[j,iel]
        isboundary(nb) && continue
        if marked[j,iel] || marked[nb.face,nb.el]
            marked[j,iel] = marked[nb.face,nb.el] = true
        end
    end
    marked
end

# Close the marks into a set that the quad patterns can refine conformingly
# (as in qmshrefine): an element with 3 marked faces gets all 4 marked, and
# otherwise the first element with 1 marked face also gets the next face in
# the cycle 1 → 3 → 2 → 4 → 1 marked. Repeat until neither case is left.
function _close_marks!(marked, m::HighOrderMesh{2,Block{2}})
    kmap = (3, 4, 2, 1)
    while true
        _share_marks!(marked, m)
        nmarked = vec(sum(marked, dims=1))
        i3 = findall(==(3), nmarked)
        if isempty(i3)
            i1 = findfirst(==(1), nmarked)
            isnothing(i1) && break
            marked[kmap[findfirst(marked[:,i1])], i1] = true
        else
            marked[:,i3] .= true
        end
    end
    marked
end

# Triangles with 2 marked faces get the third marked; 1 or 3 marked faces
# can be refined directly.
function _close_marks!(marked, m::HighOrderMesh{2,Simplex{2}})
    while true
        _share_marks!(marked, m)
        i2 = findall(==(2), vec(sum(marked, dims=1)))
        isempty(i2) && break
        marked[:,i2] .= true
    end
    marked
end

function _check_not_periodic(m::HighOrderMesh{D,G}) where {D,G}
    fmap = facemap(G())
    cel  = m.el[corner_nodes(m.fe), :]
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        nb = m.nb[j,iel]
        isboundary(nb) && continue
        sort(cel[fmap[:,j],iel]) == sort(cel[fmap[:,nb.face],nb.el]) ||
            error("Refinement of meshes with periodic faces is not supported")
    end
end

"""
    refine(m::HighOrderMesh{2}, marked)
    refine_with_parents(m::HighOrderMesh{2}, marked)

Isoparametric h-refinement of a 2D quad or triangle mesh. `marked` is a
`Bool` matrix of size `nfaces × nel` (the size of `m.nb`) marking the edges to
split; an edge marked on either side counts as marked. The marks are first
extended so that the result is conforming:

- Quads: an element with 3 marked edges gets all 4 marked, and an element
  with 1 marked edge also gets an adjacent edge marked. Two opposite marked
  edges split the element into 2 quads, two adjacent edges into 3 quads
  (with a new node at the element center), and 4 edges into 4 quads.
- Triangles: an element with 2 marked edges gets the third marked. One marked
  edge is bisected to the opposite vertex (2 triangles), and 3 marked edges
  split the element into 4 triangles.

The children take their node positions from the parent's polynomial
geometry, so curved elements stay curved and every node set and degree is
supported. Boundary tags are inherited from the parent faces. The first
child keeps the parent's element index and the other children are appended.
Meshes with periodic faces are not supported.

`refine_with_parents` also returns a vector `parent` with the index in `m` of
the parent of each element in the refined mesh.

```julia
m = set_degree(ex1mesh(), 3)
marked = falses(size(m.nb)); marked[:,1] .= true
m1 = refine(m, marked)
```
"""
refine(m::HighOrderMesh{2}, marked::AbstractMatrix{Bool}) = first(refine_with_parents(m, marked))

function refine_with_parents(m::HighOrderMesh{2,G,T}, marked::AbstractMatrix{Bool}) where {G,T}
    size(marked) == size(m.nb) ||
        throw(DimensionMismatch("marked has size $(size(marked)); expected $(size(m.nb))"))
    _check_not_periodic(m)
    marked = _close_marks!(Matrix{Bool}(marked), m)
    rules  = _refinement_rules(m.fe)

    ns, nel0 = size(m.el)
    nf       = size(m.nb, 1)
    nchildren = [ any(c) ? length(rules[c].N) : 1 for c in eachcol(marked) ]
    ntot     = sum(nchildren)

    xdg        = dg_nodes(m)
    newxdg     = Array{T}(undef, ns, ntot, 2)
    parent     = Vector{Int}(undef, ntot)
    parentface = Matrix{Int}(undef, nf, ntot)
    inew = nel0
    for iel in 1:nel0
        if !any(marked[:,iel])
            newxdg[:,iel,:]   = xdg[:,iel,:]
            parent[iel]       = iel
            parentface[:,iel] = 1:nf
            continue
        end
        rule = rules[marked[:,iel]]
        for (ic, N) in enumerate(rule.N)
            ichild = ic == 1 ? iel : (inew += 1)
            newxdg[:,ichild,:]   = N * xdg[:,iel,:]
            parent[ichild]       = iel
            parentface[:,ichild] = rule.faces[:,ic]
        end
    end

    eldg  = reshape(1:ns*ntot, ns, ntot)
    x, el = unique_mesh_nodes(reshape(newxdg, ns*ntot, 2), eldg)
    nb    = el2nb(el[corner_nodes(m.fe),:], G())
    for ichild in 1:ntot, j in 1:nf
        isboundary(nb[j,ichild]) || continue
        pf = parentface[j,ichild]
        (pf > 0 && isboundary(m.nb[pf,parent[ichild]])) ||
            error("Refined boundary face does not lie on a parent boundary face")
        nb[j,ichild] = m.nb[pf,parent[ichild]]
    end
    HighOrderMesh{2,G,T}(m.fe, x, el, nb), parent
end

"""
    uniref(m::HighOrderMesh{2}, nref=1)

Uniformly refine the 2D mesh `m` `nref` times, splitting every quad or
triangle into 4. Curved elements are refined isoparametrically, see
[`refine`](@ref).
"""
function uniref(m::HighOrderMesh{2}, nref::Integer=1)
    m = deepcopy(m)
    for _ = 1:nref
        m = refine(m, trues(size(m.nb)))
    end
    m
end

"""
    bndlayer_refine(m::HighOrderMesh{2,Block{2}}, bnds, nlayers=1)
    bndlayer_refine_with_elements(m::HighOrderMesh{2,Block{2}}, bnds, nlayers=1)

Anisotropic refinement of a quad mesh towards the boundaries with tags in
`bnds` (an integer or a collection). Each of the `nlayers` passes splits every
edge with exactly one vertex on those boundaries, so the elements along the
boundary are halved in the wall-normal direction. See [`refine`](@ref).

`bndlayer_refine_with_elements` also returns the sorted indices of the
elements in the boundary layer: those obtained by splitting an element that
had an edge marked in one of the passes.

```julia
m = rungmsh2msh("naca.geo"; porder=3)
m = bndlayer_refine(m, 1, 3)
```
"""
bndlayer_refine(m::HighOrderMesh{2,Block{2}}, bnds, nlayers::Integer=1) =
    first(bndlayer_refine_with_elements(m, bnds, nlayers))

function bndlayer_refine_with_elements(m::HighOrderMesh{2,Block{2}}, bnds, nlayers::Integer=1)
    inlayer = falses(nel(m))
    for _ = 1:nlayers
        marked = _bndlayer_marks(m, bnds)
        inlayer .|= vec(any(marked, dims=1))
        m, parent = refine_with_parents(m, marked)
        inlayer = inlayer[parent]
    end
    m, findall(inlayer)
end

# Mark the edges with exactly one vertex on the boundaries `bnds`.
function _bndlayer_marks(m::HighOrderMesh{2,Block{2}}, bnds)
    fmap  = facemap(Block{2}())
    cel   = m.el[corner_nodes(m.fe), :]
    onbnd = falses(nnodes(m))
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        nb = m.nb[j,iel]
        isboundary(nb) && bndtag(nb) ∈ bnds && (onbnd[cel[fmap[:,j],iel]] .= true)
    end
    [ count(onbnd[cel[fmap[:,j],iel]]) == 1 for j in axes(fmap,2), iel in axes(cel,2) ]
end
