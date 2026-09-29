###########################################################################
## Node deduplication

# Round x to the nearest multiple of scaling*tol to eliminate floating-point
# noise before comparisons. Adding zero() converts -0.0 → 0.0.
snap(x::T, scaling=1, tol=nothing) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling * tol * round(x / scaling / tol) + zero(T)

"""
    unique_mesh_nodes(x, el; tol=nothing, output_ix=false)

Deduplicate coincident rows in the node coordinate matrix `x` and update the
element connectivity `el` accordingly. Returns `(x, el)`, or `(x, el, ix)` if
`output_ix=true`, where `ix` maps new node indices back to rows of the original `x`.
`tol` is relative to the largest coordinate magnitude in `x`. The default is
`sqrt(eps)`, lowered when needed to stay well below the smallest distance
between two nodes of one element, so that thin elements far from the origin
(boundary layers in a large domain) keep all their nodes.
"""
function unique_mesh_nodes(x, el; tol=nothing, output_ix=false)
    scaling = maximum(abs.(x))
    tol = something(tol, _node_tol(x, el, scaling))
    xx  = snap.(x, scaling, tol)  # snap to eliminate floating-point noise
    xxx = unique(eachrow(xx))
    ix  = Int.(indexin(xxx, eachrow(xx)))  # unique row → original row
    jx  = Int.(indexin(eachrow(xx), xxx))  # original row → unique row
    x   = x[ix,:]
    el  = jx[el]
    return output_ix ? (x, el, ix) : (x, el)
end

# Default snapping tolerance of unique_mesh_nodes: sqrt(eps), or 1/1000 of the
# smallest max-norm distance between two distinct nodes of one element
# (relative to scaling) if that is smaller, but not below 100 eps.
function _node_tol(x, el, scaling)
    T = float(eltype(x))
    dmin = T(Inf)
    for iel in axes(el,2), i in axes(el,1), j in 1:i-1
        a, b = el[i,iel], el[j,iel]
        d = maximum(abs(x[a,k] - x[b,k]) for k in axes(x,2))
        d > 0 && (dmin = min(dmin, d))
    end
    scaling > 0 && isfinite(dmin) || return sqrt(eps(T))
    clamp(dmin / scaling / 1000, 100eps(T), sqrt(eps(T)))
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
        nb = m.nb[j,iel]
        if isboundary(nb) && (isnothing(bndnbrs) || bndtag(nb) ∈ bndnbrs)
            append!(nodes, m.el[f2n[:,j], iel])
        end
    end
    unique(nodes)
end

"""
    boundary_distance(m::HighOrderMesh, bndnbrs=nothing; nsub=32)

Distance from every DG node of `m` to the nearest boundary face, as an
`ns × nel` matrix in the node layout of [`dg_nodes`](@ref). If `bndnbrs` is
given (an integer or collection of integers), only faces with those boundary
numbers count; otherwise all boundary faces do. Typical uses are the wall
distance of turbulence models and sizing functions near a boundary.

Each curved face is sampled on a uniform grid of its reference element with
`nsub` intervals per edge and replaced by the segments (2D) or triangles (3D)
between the samples. The distance is exact for straight faces, and otherwise
off by about `κ h^2 / 8` for a face of size `h` and curvature `κ`. Faces whose
bounding box is farther away than the closest face found so far are skipped,
but there is no spatial search structure: the cost grows like the number of
nodes times the number of boundary faces.

```julia
msh = mshairfoil(:naca0012; aoa=0)
d   = boundary_distance(msh, 1; nsub=128)   # distance to the airfoil
```
"""
function boundary_distance(m::HighOrderMesh{D,G,T}, bndnbrs=nothing; nsub::Integer=32) where {D,G,T}
    D in (2, 3) || error("boundary_distance supports 2D and 3D meshes")
    f2n = mkface2nodes(m)
    fe  = subelement(m.fe, D-1)
    ξ, simplices = _face_samples(elgeom(fe), nsub)
    N   = shapefcns(fe, ξ)
    faces = Vector{SVector{D,T}}[]
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        nb = m.nb[j,iel]
        if isboundary(nb) && (isnothing(bndnbrs) || bndtag(nb) ∈ bndnbrs)
            xs = N * m.x[m.el[f2n[:,j],iel],:]
            push!(faces, [ SVector{D,T}(r) for r in eachrow(xs) ])
        end
    end
    isempty(faces) && error("No boundary faces with the numbers $bndnbrs")
    lo = [ reduce((a,b) -> min.(a,b), X) for X in faces ]
    hi = [ reduce((a,b) -> max.(a,b), X) for X in faces ]
    boxdist(p, f) = norm(max.(lo[f] .- p, p .- hi[f], zero(T)))
    facedist(p, f) = minimum(s -> _simplex_distance(p, map(k -> faces[f][k], s)...), simplices)

    x = dg_nodes(m)
    d = zeros(T, size(x,1), size(x,2))
    for iel in axes(x,2), i in axes(x,1)
        p = SVector{D,T}(ntuple(k -> x[i,iel,k], D))
        dmin = facedist(p, argmin(f -> boxdist(p, f), eachindex(faces)))
        for f in eachindex(faces)
            boxdist(p, f) < dmin && (dmin = min(dmin, facedist(p, f)))
        end
        d[i,iel] = dmin
    end
    d
end

# Uniform sample points ξ (rows) on a face reference element with n intervals
# per edge, in the canonical node order (first index fastest), and the segments
# or triangles between them as tuples of row indices.
_face_samples(::Block{1}, n) = collect(reshape((0:n) ./ n, :, 1)), [ (i, i+1) for i in 1:n ]

function _face_samples(::Block{2}, n)
    ix(i, j) = 1 + i + (n+1)*j
    ξ = [ ij[c] / n for ij in [ (i, j) for j in 0:n for i in 0:n ], c in 1:2 ]
    tris = NTuple{3,Int}[]
    for j in 0:n-1, i in 0:n-1
        push!(tris, (ix(i,j), ix(i+1,j), ix(i+1,j+1)), (ix(i,j), ix(i+1,j+1), ix(i,j+1)))
    end
    ξ, tris
end

function _face_samples(::Simplex{2}, n)
    ijs = [ (i, j) for j in 0:n for i in 0:n-j ]
    ix  = Dict(ij => k for (k, ij) in enumerate(ijs))
    ξ   = [ ij[c] / n for ij in ijs, c in 1:2 ]
    tris = NTuple{3,Int}[]
    for (i, j) in ijs
        i + j <= n-1 && push!(tris, (ix[(i,j)], ix[(i+1,j)], ix[(i,j+1)]))
        i + j <= n-2 && push!(tris, (ix[(i+1,j)], ix[(i+1,j+1)], ix[(i,j+1)]))
    end
    ξ, tris
end

# Distance from p to the segment ab
function _simplex_distance(p, a, b)
    ab = b - a
    l2 = dot(ab, ab)
    t  = l2 > 0 ? clamp(dot(p - a, ab) / l2, 0, 1) : zero(l2)
    norm(p - a - t * ab)
end

# Distance from p to the triangle abc, with the closest point found from the
# Voronoi regions of the vertices and edges (Ericson, Real-Time Collision
# Detection, 2005)
function _simplex_distance(p, a, b, c)
    ab, ac, ap = b - a, c - a, p - a
    d1, d2 = dot(ab, ap), dot(ac, ap)
    d1 <= 0 && d2 <= 0 && return norm(ap)
    bp = p - b
    d3, d4 = dot(ab, bp), dot(ac, bp)
    d3 >= 0 && d4 <= d3 && return norm(bp)
    vc = d1*d4 - d3*d2
    vc <= 0 && d1 >= 0 && d3 <= 0 && return _simplex_distance(p, a, b)
    cp = p - c
    d5, d6 = dot(ab, cp), dot(ac, cp)
    d6 >= 0 && d5 <= d6 && return norm(cp)
    vb = d5*d2 - d1*d6
    vb <= 0 && d2 >= 0 && d6 <= 0 && return _simplex_distance(p, a, c)
    va = d3*d6 - d5*d4
    va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0 && return _simplex_distance(p, b, c)
    den = va + vb + vc
    den > 0 || return min(_simplex_distance(p, a, b), _simplex_distance(p, a, c),
                          _simplex_distance(p, b, c))    # degenerate triangle
    norm(ap - (vb/den) * ab - (vc/den) * ac)
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
        isboundary(m.nb[j,iel]) || continue  # interior face
        facex  = m.x[m.el[f2n[:,j],iel],:]
        onbnd  = hcat([ snap.(bndexpr(cx)) .== 0 for cx in eachrow(facex) ]...)
        bndnbr = findfirst(all(onbnd, dims=2)[:])
        isnothing(bndnbr) && error("No boundary expression matching boundary face")
        m.nb[j,iel] = Neighbor(-bndnbr, 0, 0)
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
"""
function set_bnd_periodic!(m::HighOrderMesh{D,G,T}, bnds, dir) where {D,G,T}
    f2n  = mkface2nodes(m)
    fmap = facemap(G())
    corner_el = m.el[corner_nodes(m.fe), :]
    match_coords = (1:D) .≠ dir  # compare all coords except the periodic direction

    function periodic_face_permutation(x1, x2)
        D <= 2 && return Int16(1)
        sx1 = snap.(x1[:,match_coords])
        sx2 = snap.(x2[:,match_coords])
        perm = findfirst(i -> sx2[i,:] == sx1[1,:], axes(sx2,1))
        isnothing(perm) && error("Could not determine periodic face permutation")
        Int16(perm)
    end

    dd = Dict{Matrix{T}, NTuple{2,Int}}()
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        bndtag(m.nb[j,iel]) ∈ bnds || continue
        facex = m.x[m.el[f2n[:,j],iel],:]
        key   = sortslices(snap.(facex[:,match_coords]), dims=1)
        if haskey(dd, key)
            iel0, j0 = pop!(dd, key)
            fcx  = m.x[corner_el[fmap[:,j],iel],:]
            fcx0 = m.x[corner_el[fmap[:,j0],iel0],:]
            m.nb[j,iel]   = Neighbor(iel0, j0, periodic_face_permutation(fcx, fcx0))
            m.nb[j0,iel0] = Neighbor(iel,  j,  periodic_face_permutation(fcx0, fcx))
        else
            dd[key] = (iel, j)
        end
    end
end
