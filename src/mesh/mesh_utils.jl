###########################################################################
## Uniform h-refinement

"""
    uniref(m::HighOrderMesh, nref=1)

Uniformly refine `m` by splitting each element into 4 (2D triangles/quads)
subsets. Currently implemented for 2D linear (`p=1`) meshes only.
Apply `nref` times by passing an integer.
"""
function uniref(m::HighOrderMesh, nref)
    m = deepcopy(m)
    for _ = 1:nref
        m = uniref(m)
    end
    m
end

function uniref(m::HighOrderMesh{2,G,1,T}) where {G,T}
    emap   = edgemap(G())
    # Collect and deduplicate all element edges (sorted vertex pairs)
    eledges = sort.([ cel[edge] for edge = eachcol(emap), cel = eachcol(m.el) ])
    edges   = unique(eledges)
    eidx    = indexin(eledges, edges)  # edge index for each (face, element) pair

    nx     = size(m.x, 1)
    nel    = size(m.el, 2)
    nedges = length(edges)

    # New nodes: edge midpoints, and (for quads) element centroids
    pmid   = [ sum(m.x[edge,d])/2 for edge in edges, d = 1:2 ]
    mapmid = eidx .+ nx  # global indices of edge midpoints

    mkaddx(::Simplex) = zeros(0, 2)
    mkaddx(::Block)   = [ sum(m.x[el,d])/4 for el = eachcol(m.el), d = 1:2 ]
    newx = vcat(m.x, pmid, mkaddx(G()))

    # Build sub-element connectivity (geometry-specific)
    function mkels(::Simplex)
        prev = [3,1,2]
        next = [2,3,1]
        # 3 corner sub-triangles + 1 central sub-triangle per element
        ts  = [ [m.el[i,it], mapmid[prev[i],it], mapmid[next[i],it]] for i = 1:3, it = 1:nel ]
        ts2 = [ mapmid[:,it] for it = 1:nel ]
        vcat(ts[:], ts2[:])
    end

    function mkels(::Block)
        mapc = nx + nedges .+ (1:nel)  # centroid node indices
        prev = [4,1,2,3]
        [ [m.el[emap[1,i],:] mapmid[i,:] mapmid[prev[i],:] mapc]' for i = 1:4 ]
    end

    newel = hcat(mkels(G())...)
    return HighOrderMesh(newx, newel)
end

###########################################################################
## Node deduplication

# Round x to the nearest multiple of scaling*tol to eliminate floating-point
# noise before comparisons. Adding zero() converts -0.0 → 0.0.
snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling * tol * round(x / scaling / tol) + zero(T)

"""
    unique_mesh_nodes(x, el; output_ix=false)

Deduplicate coincident rows in the node coordinate matrix `x` and update the
element connectivity `el` accordingly. Returns `(x, el)`, or `(x, el, ix)` if
`output_ix=true`, where `ix` maps new node indices back to rows of the original `x`.
"""
function unique_mesh_nodes(x, el; output_ix=false)
    xx  = snap.(x, maximum(abs.(x)))  # snap to eliminate floating-point noise
    xxx = unique(eachrow(xx))
    ix  = Int.(indexin(xxx, eachrow(xx)))  # unique row → original row
    jx  = Int.(indexin(eachrow(xx), xxx))  # original row → unique row
    x   = x[ix,:]
    el  = jx[el]
    return output_ix ? (x, el, ix) : (x, el)
end

###########################################################################
## Boundary conditions

"""
    boundary_nodes(m::HighOrderMesh, bndnbrs=nothing)

Return the global node indices on boundary faces. If `bndnbrs` is given
(an integer or collection of integers), only faces with those boundary
numbers are included; otherwise all boundary faces are returned.
"""
function boundary_nodes(m::HighOrderMesh, bndnbrs=nothing)
    f2n  = mkface2nodes(m)
    nf, nel = size(m.nb)
    nodes = Int64[]
    for iel in 1:nel, j in 1:nf
        jel, _, _ = m.nb[j,iel]
        if jel < 1 && (isnothing(bndnbrs) || -jel ∈ bndnbrs)
            append!(nodes, m.el[f2n[:,j], iel])
        end
    end
    unique(nodes)
end

"""
    set_bnd_numbers!(m::HighOrderMesh, bndexpr)

Label each boundary face in `m.nb` with a boundary region number.
`bndexpr` is a function `x -> [expr1(x), expr2(x), ...]` where `expri(x) == 0`
for all nodes on boundary region `i`. Errors if a face matches no expression.
"""
function set_bnd_numbers!(m::HighOrderMesh, bndexpr)
    f2n     = mkface2nodes(m)
    nf, nel = size(m.nb)
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        m.nb[j,iel][1] >= 1 && continue  # interior face
        facex  = m.x[m.el[f2n[:,j],iel],:]
        onbnd  = hcat([ snap.(bndexpr(cx)) .== 0 for cx in eachrow(facex) ]...)
        bndnbr = findfirst(all(onbnd, dims=2)[:])
        isnothing(bndnbr) && error("No boundary expression matching boundary face")
        m.nb[j,iel] = (-bndnbr, 0, 0)
    end
end

"""
    set_bnd_periodic!(m::HighOrderMesh, bnds, dir)

Connect the boundary faces numbered in `bnds` periodically along coordinate
direction `dir`. Matching is done by comparing the face node coordinates in
all directions except `dir`.

```julia
msh = mshsquare(5)
set_bnd_periodic!(msh, (1,2), 1)   # periodic left/right (x)
set_bnd_periodic!(msh, (3,4), 2)   # periodic bottom/top (y)
```

TODO: neighbor face permutation (orientation) is not yet tracked.
"""
function set_bnd_periodic!(m::HighOrderMesh{D,G,P,T}, bnds, dir) where {D,G,P,T}
    f2n  = mkface2nodes(m)
    match_coords = (1:D) .≠ dir  # compare all coords except the periodic direction
    dd = Dict{Matrix{T}, NTuple{2,Int}}()
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        -m.nb[j,iel][1] ∈ bnds || continue
        facex = m.x[m.el[f2n[:,j],iel],:]
        key   = sortslices(snap.(facex[:,match_coords]), dims=1)
        if haskey(dd, key)
            iel0, j0 = pop!(dd, key)
            m.nb[j,iel]   = (iel0, j0, 0)
            m.nb[j0,iel0] = (iel,  j,  0)
        else
            dd[key] = (iel, j)
        end
    end
end
