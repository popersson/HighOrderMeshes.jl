using HDF5

######################
## File I/O

filevars = ("s1","x","el","nb")

function savemesh(fname, m::HighOrderMesh)
    h5open(fname, "w") do fid
        dat = (m.fe.ref_nodes[1][:], m.x, m.el, m.nb)
        for i in eachindex(filevars)
            fid[filevars[i], compress=9] = dat[i]
        end
    end
end

function loadmesh(fname)
    s1,x,el,nb = h5open(fname, "r") do fid
        read.([fid], filevars)
    end
    nb = NeighborData.(nb)
    D = size(x,2)
    P = length(s1) - 1
    G = find_elgeom(D, P, size(el,1))
    fe = FiniteElement(G, s1)
    m = HighOrderMesh(fe, x, el, nb)
end

######################


"""
    write_matrix(f, x)

TBW
"""
function write_matrix(f, x)
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            print(f, x[i,j])
            if j<size(x,2)
                print(f, ' ')
            end
        end
        println(f)
    end
end

"""
    read_matrix!(f, x)

TBW
"""
function read_matrix!(f, x)
    for i = 1:size(x,1)
        line = readline(f)
        x[i,:] = parse.(eltype(x), split(line, ' ')[1:size(x,2)])
    end
end

"""
    savemeshtxt(fname, m::HighOrderMesh)

TBW
"""
function savemeshtxt(fname, m::HighOrderMesh)
    nx = size(m.x,1)
    nel = size(m.el,2)

    f = open(fname, "w")
    println(f, "$nx $nel $(size(m.el,1))  # nbr_nodes nbr_elems nbr_nodes_per_elem")
    write_matrix(f, m.x)
    write_matrix(f, m.el')
    close(f)
end

"""
    loadmeshtxt(fname)

TBW
"""
function loadmeshtxt(fname)
    f = open(fname, "r")
    line = split(readline(f), " ")
    nx,nel,ne = parse.(Int64, line[1:3])
    name = line[3]
    x = zeros(Float64, nx, 2)
    read_matrix!(f, x)
    el = zeros(Int64, nel, ne)
    read_matrix!(f, el)
    close(f)
    return HighOrderMesh(x,el')
end

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
    mkface2nodes(eg, sface::AbstractArray{T}, svol::AbstractArray{T}) where T
    mkface2nodes(fe::FiniteElement)

TBW
"""
function mkface2nodes(eg::ElementGeometry{D}, sface::AbstractArray{T}, svol::AbstractArray{T}) where {D,T}
    fmap = facemap(eg)
    basis_face = snap.(eval_shapefcns(subgeom(eg,D-1), sface))
    basis_vol = snap.(eval_shapefcns(eg, svol))
    f2n = fill(0, size(basis_face,1), size(fmap,2))
    for (ii,ic) in enumerate(eachcol(fmap))
        ix = indexin(eachrow(basis_face), eachrow(basis_vol[:,ic]))
        f2n[:,ii] .= ix
    end
    f2n
end

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

mkface2nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T} =
    mkface2nodes(G(), ref_nodes(fe, D-1), ref_nodes(fe, D))

mkface2nodes(m::HighOrderMesh) = mkface2nodes(m.fe)

function mkldgswitch(eg::Block{D}, nb) where {D}
    nf,nel = size(nb)
    sw = fill(-1, nf, nel)

    opposite_face(i) = i - iseven(i) + isodd(i)

    while true
        idx = findfirst(isequal(-1), sw)
        isnothing(idx) && break
        j,iel = Tuple(idx)

        iel0 = iel
        j0 = j
        dir = 1
        while true
            sw[j,iel] = dir
            jel,k,_ = nb[j,iel]
            if jel > 0
                sw[k,jel] = 1 - dir
                iel = jel
                j = opposite_face(k)
                if sw[j,iel] != -1
                    if dir == 1
                        dir = 0
                        iel = iel0
                        j = opposite_face(j0)
                    else
                        break
                    end
                end
            elseif dir == 1
                dir = 0
                iel = iel0
                j = opposite_face(j0)
            else
                break
            end
        end
    end
    return sw
end

mkldgswitch(m::HighOrderMesh) = mkldgswitch(elgeom(m), m.nb)

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

function align_with_ldgswitch!(m::HighOrderMesh{2,Block{2},P}, sw=nothing) where P
    if isnothing(sw)
        sw = mkldgswitch(m)
    end

    switch_cases = [ [1,0,1,0],
                     [1,0,0,1],
                     [0,1,0,1],
                     [0,1,1,0] ]

    mapcase = [ findfirst([csw] .== switch_cases) for csw in eachcol(sw) ]

    nds = collect(reshape(1:(P+1)^2, P+1, P+1))
    ndmaps = [nds]
    for i = 2:4
        push!(ndmaps, rotl90(ndmaps[end]))
    end
    ndmaps = reshape.(ndmaps, :)
    
    fcmaps = [ [1,2,3,4],
               [4,3,1,2],
               [2,1,4,3],
               [3,4,2,1] ]
    ifcmaps = [ invperm(fc) for fc in fcmaps ]
    
    for iel = 1:size(m.el,2)
        cmap = mapcase[iel]
        m.el[:,iel] .= m.el[:,iel][ndmaps[cmap]]
        cnb = m.nb[:,iel]
        for j = 1:4
            if cnb[j][2] > 0
                jel = cnb[j][1]
                cnb[j] = (jel, ifcmaps[mapcase[jel]][cnb[j][2]], 0 )
            end
        end
        m.nb[:,iel] = cnb[fcmaps[cmap]]
    end
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
            isnothing(bndnbr) && throw("No boundary expression matching boundary face")
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

