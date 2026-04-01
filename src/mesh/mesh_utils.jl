"""
    uniref(m::HighOrderMesh, nref)
    uniref(m::HighOrderMesh{2,G,1,T}) where {G,T}

TBW
"""
function uniref(m::HighOrderMesh, nref)
    m = deepcopy(m)
    for i = 1:nref
        m = uniref(m)
    end
    m
end

function uniref(m::HighOrderMesh{2,G,1,T}) where {G,T}
    emap = edgemap(G())
    eledges = sort.([ cel[edge] for edge = eachcol(emap), cel = eachcol(m.el) ])
    edges = unique(eledges)
    map = indexin(eledges, edges)

    nx = size(m.x,1)
    nel = size(m.el,2)
    nedges = length(edges)

    pmid = [ sum(m.x[edge,dim])/2 for edge in edges, dim = 1:2 ]
    mapmid = map .+ nx

    mkaddx(::Simplex) = zeros(0,2)
    mkaddx(::Block) = [ sum(m.x[el,dim])/4 for el = eachcol(m.el), dim = 1:2 ]
    newx = vcat(m.x, pmid, mkaddx(G()))

    function mkels(::Simplex)
        prev = [3,1,2]
        next = [2,3,1]
        ts = [ [m.el[i,it],mapmid[prev[i],it],mapmid[next[i],it]] for i = 1:3, it = 1:nel ]
        ts2 = [ mapmid[:,it] for it = 1:nel ]
        els = vcat(ts[:], ts2[:])
    end

    function mkels(::Block)
        mapc = nx + nedges .+ (1:nel)
        prev = [4,1,2,3]
        els = [ [m.el[emap[1,i],:] mapmid[i,:] mapmid[prev[i],:] mapc]' for i = 1:4 ]
    end

    newel = hcat(mkels(G())...)

    return HighOrderMesh(newx,newel)
end

snap(x::T, scaling=1) where {T <: Real} = x
snap(x::T, scaling=1, tol=sqrt(eps(T))) where {T <: AbstractFloat} =
    scaling*tol*round(x/scaling/tol) + zero(T)  # Adding zero to uniquify -0.0 and 0.0

"""
    unique_mesh_nodes(x, el)

TBW
"""
function unique_mesh_nodes(x, el; output_ix=false)
    xx = snap.(x, maximum(abs.(x)))
    xxx = unique(eachrow(xx))
    ix = Int.(indexin(xxx, eachrow(xx)))
    jx = Int.(indexin(eachrow(xx), xxx))
    x = x[ix,:]
    el = jx[el]
    return output_ix ? (x,el,ix) : (x,el)
end

function boundary_nodes(m::HighOrderMesh, bndnbrs=nothing)
    f2n = mkface2nodes(m)
    nf,nel = size(m.nb)
    
    nodes = Int64[]
    for iel in 1:nel
        for j in 1:nf
            jel,k,_ = m.nb[j,iel]
            if jel < 1 && (isnothing(bndnbrs) || -jel ∈ bndnbrs)
                append!(nodes, m.el[f2n[:,j],iel])
            end
        end
    end
    unique(nodes)
end

function set_bnd_numbers!(m::HighOrderMesh, bndexpr)
    nf,nel = size(m.nb)
    f2n = mkface2nodes(m)

    scaling = maximum(abs.(m.x))
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        if m.nb[j,iel][1] < 1
            facex = m.x[m.el[f2n[:,j],iel],:]
            onbnd = hcat([ snap.(bndexpr(cx)) .== 0 for cx in eachrow(facex) ]...)
            bndnbr = findfirst(all(onbnd,dims=2)[:])
            isnothing(bndnbr) && error("No boundary expression matching boundary face")
            m.nb[j,iel] = (-bndnbr,0,0)
        end
    end
end

"""
    set_bnd_periodic!(m::HighOrderMesh{D,G,P,T}, bnds, dir) where {D,G,P,T}

Example:
```julia
    msh = mshsquare(5)
    set_bnd_periodic!(msh, (1,2), 1)    # Periodic left/right (x-direction)
    set_bnd_periodic!(msh, (3,4), 2)    # Periodic bottom/top (y-direction)
```
"""
function set_bnd_periodic!(m::HighOrderMesh{D,G,P,T}, bnds, dir) where {D,G,P,T}
    # TODO: Neighbor face permutation
    nf,nel = size(m.nb)
    f2n = mkface2nodes(m)

    match_coords = (1:D) .≠ dir
    scaling = maximum(abs.(m.x))
    dd = Dict{Matrix{T}, NTuple{2,Int}}()
    for iel in axes(m.nb,2), j in axes(m.nb,1)
        if -m.nb[j,iel][1] ∈ bnds
            facex = m.x[m.el[f2n[:,j],iel],:]
            key = snap.(facex[:,match_coords])
            key = sortslices(key, dims=1)
            if haskey(dd, key)
                iel0,j0 = pop!(dd, key)
                m.nb[j,iel] = (iel0,j0,0)
                m.nb[j0,iel0] = (iel,j,0)
            else
                dd[key] = (iel,j)
            end
        end
    end
end
