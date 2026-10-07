###########################################################################
## Curved boundaries from a signed distance function

"""
    project_points(x, fd; maxiter=10)

Project the points in the rows of `x` (`npoints × D`) onto the zero level set
of the signed distance function `fd` and return them as a new matrix. `fd(p)`
takes a point `p` as an `SVector{D}` (use `p[1]`, `p[2]`, `norm(p)`, ...) and
is negative inside the domain, as in DistMesh.

Each point takes Newton steps `p ← p - fd(p) ∇fd(p) / |∇fd(p)|²`, with the
gradient by central differences, until `|fd(p)|` is at round-off level or
after `maxiter` steps. For an exact signed distance function this is the
closest point on the boundary; for any other function with the same zero
level set it is a nearby point on it.

```julia
fd(p) = sqrt(sum(p.^2)) - 1                   # unit circle
x = project_points([2.0 0.0; 0.5 0.5], fd)    # [1 0; √½ √½]
```
"""
function project_points(x::AbstractMatrix, fd; maxiter::Integer=10)
    T = float(eltype(x))
    isempty(x) && return Matrix{T}(x)
    _project_points(x, fd, _point_scale(x), maxiter, Val(size(x,2)))
end

# Length scale of a point set for finite difference steps and tolerances
function _point_scale(x)
    T = float(eltype(x))
    extent = maximum(maximum(c) - minimum(c) for c in eachcol(x))
    s = max(extent, maximum(abs, x))
    s > 0 ? T(s) : one(T)
end

function _project_points(x, fd, scale::T, maxiter, ::Val{D}) where {T,D}
    h   = cbrt(eps(T)) * scale
    tol = 100 * eps(T) * scale
    y   = Matrix{T}(undef, size(x))
    for i in axes(x,1)
        p = SVector{D,T}(ntuple(k -> x[i,k], D))
        for _ in 1:maxiter
            d = fd(p)
            abs(d) <= tol && break
            g = SVector{D,T}(ntuple(k -> (fd(p + h*_unit(Val(D), k, T)) -
                                          fd(p - h*_unit(Val(D), k, T))) / 2h, D))
            g2 = dot(g, g)
            g2 > 0 || break
            p = p - T(d / g2) * g
        end
        y[i,:] = p
    end
    y
end

_unit(::Val{D}, k, T) where {D} = SVector{D,T}(ntuple(j -> j == k ? one(T) : zero(T), D))

"""
    interp_arclength(x, s)

Points at the fractions `s` (in `[0,1]`) of the arclength of the polyline
through the rows of `x` (`npoints × D`), as a `length(s) × D` matrix. The
points are interpolated linearly between the rows of `x`; `s = 0` and
`s = 1` give the first and last rows exactly.

```julia
x = [0.0 0.0; 1.0 0.0; 1.0 3.0]                # polyline of length 4
interp_arclength(x, [0, 0.25, 0.5, 1])         # [0 0; 1 0; 1 1; 1 3]
```
"""
function interp_arclength(x::AbstractMatrix, s::AbstractVector)
    T = float(promote_type(eltype(x), eltype(s)))
    n = size(x,1)
    n >= 1 || throw(ArgumentError("the polyline needs at least one point"))
    L = zeros(T, n)
    for i in 2:n
        L[i] = L[i-1] + norm(view(x,i,:) - view(x,i-1,:))
    end
    y = Matrix{T}(undef, length(s), size(x,2))
    for (k, sk) in enumerate(s)
        if n == 1 || L[end] == 0
            y[k,:] = view(x,1,:)
            continue
        end
        target = sk * L[end]
        i = clamp(searchsortedlast(L, target), 1, n-1)
        len = L[i+1] - L[i]
        α = len > 0 ? (target - L[i]) / len : zero(T)
        y[k,:] = (1 - α) * view(x,i,:) + α * view(x,i+1,:)
    end
    y
end

"""
    curve_boundary(m::HighOrderMesh{2}, fd, bndnbrs=nothing; maxiter=20)

Return a new mesh with the boundary faces of `m` moved onto the zero level
set of the signed distance function `fd` (negative inside, called as
`fd(p)` with an `SVector` `p`, as in DistMesh). If `bndnbrs` is given (an
integer or collection of integers), only faces with those boundary numbers
move. This is the usual way to turn a straight-sided mesh, for example from
DistMesh, into a curved high-order mesh: raise the degree with
[`set_degree`](@ref) first, then curve the boundary.

1. The element vertices on the faces are projected onto the curve with
   [`project_points`](@ref).
2. The other nodes of each face are placed on the curve at the arclength
   fractions of the reference face nodes, by alternating
   [`interp_arclength`](@ref) and `project_points` until they settle (at
   most `maxiter` times), with the vertices held fixed.
3. The displacements of the face nodes are extended into the elements by
   blending: transfinite (Coons) interpolation of the four edges for quads,
   and the polynomial edge blending of Szabó and Babuška for triangles.
   Elements that only touch a moved vertex follow it linearly, so the mesh
   stays continuous, and the new node positions are polynomial in the
   reference coordinates. Straight faces that already lie on the zero level
   set do not move.

Works for any degree and any reference node set of quads and triangles.
Faces may not be too coarse for the curvature of the boundary: the blending
does not check for inverted elements.

```julia
using DistMesh                                 # optional, for distmesh2d
fd(p) = sqrt(sum(p.^2)) - 1                    # unit circle
msh = HighOrderMesh(distmesh2d(fd, huniform, 0.2, ((-1,-1), (1,1))))
msh = curve_boundary(set_degree(msh, 4), fd)   # degree 4, curved boundary
```
"""
function curve_boundary(m::HighOrderMesh{2,G,T}, fd, bndnbrs=nothing; maxiter::Integer=20) where {G,T}
    fe    = m.fe
    f2n   = mkface2nodes(m)
    fmap  = facemap(G())
    cn    = corner_nodes(fe)
    iscurved = [ isboundary(nb) && (isnothing(bndnbrs) || bndtag(nb) ∈ bndnbrs) for nb in m.nb ]
    scale = _point_scale(m.x)
    proj(y) = _project_points(y, fd, scale, maxiter, Val(2))

    # 1. Vertices of the curved faces
    x   = copy(m.x)
    vix = unique(m.el[cn[fmap[k,j]], iel] for k in 1:2, (j, iel) in Tuple.(findall(iscurved)))
    x[vix,:] = proj(x[vix,:])

    # 2. Face nodes on the curve, at the arclength fractions of the reference face nodes
    s    = vec(ref_nodes(fe, 1))
    perm = sortperm(s)
    for I in findall(iscurved)
        j, iel = Tuple(I)
        ix = m.el[f2n[perm,j], iel]
        X  = x[ix,:]
        for _ in 1:maxiter
            X0 = X
            X  = proj(interp_arclength(X, s[perm]))
            X[[1,end],:] = X0[[1,end],:]
            maximum(abs, X - X0) <= 1000 * eps(T) * scale && break
        end
        x[ix,:] = X
    end

    # 3. Blend the face displacements into the elements
    δ = x - m.x
    N, Nf, B = _boundary_blending(fe)
    for iel in axes(m.el,2)
        nodes = view(m.el, :, iel)
        any(!iszero, view(δ, nodes, :)) || continue
        δv  = δ[nodes[cn], :]
        δel = N * δv
        for j in axes(m.nb,1)
            iscurved[j,iel] || continue
            φ = δ[nodes[f2n[:,j]], :] - Nf * δv[fmap[:,j], :]   # face displacement minus its linear part
            δel += B[j] * φ
        end
        x[nodes,:] = m.x[nodes,:] + δel
    end
    HighOrderMesh{2,G,T}(fe, x, copy(m.el), copy(m.nb))
end

# Matrices for the displacement blending of curve_boundary: N (ns × nv) maps
# vertex displacements to the nodes, Nf (nfacenodes × 2) maps the displacements
# of a face's two vertices to its nodes, and B[j] (ns × nfacenodes) maps the
# nodal values of a displacement of face j that vanishes at the vertices to
# the nodes of the element.
function _boundary_blending(fe::FiniteElement{2,G,T}) where {G,T}
    ξ    = ref_nodes(fe)
    N    = shapefcns(G(), ξ)
    fe1  = subelement(fe, 1)
    Nf   = shapefcns(elgeom(fe1), ref_nodes(fe1))
    fmap = facemap(G())
    B = map(axes(fmap,2)) do j
        tw = [ _edge_blend(G(), ξ[i,:], N[i,:], fmap[1,j], fmap[2,j]) for i in axes(ξ,1) ]
        last.(tw) .* shapefcns(fe1, reshape(first.(tw), :, 1))
    end
    N, Nf, B
end

# Edge parameter t (0 at vertex a, 1 at vertex b) and blending weight of the
# edge from vertex a to vertex b at the reference point ξ with linear shape
# function values λ. Quads: Coons weight 1 - (distance to the edge).
function _edge_blend(::Block{2}, ξ, _, a, b)
    V  = vertices(Block{2}())
    e  = V[b,:] - V[a,:]
    t  = dot(ξ - V[a,:], e) / dot(e, e)
    t, 1 - norm(ξ - V[a,:] - t * e)
end

# Triangles: Szabó-Babuška blending λa λb φ(t) / (t (1-t)) with t = (1 + λb - λa) / 2,
# a polynomial of the same degree as the edge displacement φ, which vanishes at t = 0, 1.
function _edge_blend(::Simplex{2}, _, λ, a, b)
    t   = (1 + λ[b] - λ[a]) / 2
    den = t * (1 - t)
    t, den > 0 ? λ[a] * λ[b] / den : zero(den)
end
