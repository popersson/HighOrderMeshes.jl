###########################################################################
## Node deduplication

# Round x to the nearest multiple of scaling*tol to eliminate floating-point
# noise before comparisons. Adding zero() converts -0.0 → 0.0.
snap(x::T, scaling=1, tol=nothing) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling * tol * round(x / scaling / tol) + zero(T)

"""
    unique_mesh_nodes(x, el; tol=sqrt(eps(T)), output_ix=false)

Deduplicate coincident rows in the node coordinate matrix `x` and update the
element connectivity `el` accordingly. Returns `(x, el)`, or `(x, el, ix)` if
`output_ix=true`, where `ix` maps new node indices back to rows of the original `x`.
`tol` is relative to the largest coordinate magnitude in `x`.
"""
function unique_mesh_nodes(x, el; tol=sqrt(eps(float(eltype(x)))), output_ix=false)
    xx  = snap.(x, maximum(abs.(x)), tol)  # snap to eliminate floating-point noise
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
        nb = m.nb[j,iel]
        if isboundary(nb) && (isnothing(bndnbrs) || bndtag(nb) ∈ bndnbrs)
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
