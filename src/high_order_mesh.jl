struct HighOrderMesh{D,G,P,T}
    fe::FiniteElement{D,G,P,T}
    x::Matrix{T}
    el::Matrix{Int}
    nbor::Matrix{NTuple{2,Int}}
end

"""
    HighOrderMesh(fe::FiniteElement{D,G,P,T},
              x::AbstractMatrix{T},
              el::AbstractMatrix{Int}) where {D,G,P,T}
    HighOrderMesh(x::Matrix{T}, el::AbstractMatrix{Int}) where {T}

TBW
"""
function HighOrderMesh(fe::FiniteElement{D,G,P,T},
              x::AbstractMatrix{T},
              el::AbstractMatrix{Int}) where {D,G,P,T}
    nbor = el2nbor(el[corner_nodes(fe),:], G())
    HighOrderMesh{D,G,P,T}(fe, x, el, nbor)
end

function HighOrderMesh(x::Matrix{T}, el::AbstractMatrix{Int}) where {T}
    dim, nv = size(x,2), size(el,1)
    eg = find_elgeom(dim, nv)
    nbor = el2nbor(el, eg)
    fe = FiniteElement(eg, 1, T)
    HighOrderMesh{dim,typeof(eg),1,T}(fe, x, el, nbor)
end

dg_nodes(m::HighOrderMesh) = m.x[m.el,:]

elgeom(::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = G()
dim(::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = D
porder(m::HighOrderMesh{D,G,P,T}) where {D,G,P,T} = P

function Base.show(io::IO, m::HighOrderMesh)
    print(io, "HighOrderMesh: $(dim(m))D, ")
    print(io, "$(size(m.x,1)) vertices, ")
    print(io, "$(size(m.el,2)) $(name(m.fe)) elements.")
end

function el2nbor(el, eg)
    map = facemap(eg)
    nv,nel = size(el)
    nfv = size(map,1)
    nf = size(map,2)

    nb = fill((0,0), nf, nel)
    dd = Dict{NTuple{nfv,Int}, NTuple{2,Int}}()
    sizehint!(dd, nel*nf)
    e = fill(0,nfv)
    for iel = 1:nel
        for jf = 1:nf
            e[:] = el[map[:,jf],iel]
            sort!(e)
            et = Tuple(e)
            if haskey(dd,et)
                nbel = pop!(dd,et)
                nb[jf,iel] = nbel
                nb[nbel[2],nbel[1]] = (iel,jf)
            else
                dd[et] = (iel,jf)
            end
        end
    end
    nb
end

"""
    change_ref_nodes(m::HighOrderMesh{D,G,P,T}, newfe::FiniteElement) where {D,G,P,T}

TBW
"""
function change_ref_nodes(m::HighOrderMesh{D,G,P,T}, newfe::FiniteElement) where {D,G,P,T}
    newp = porder(newfe)
    Pfe = eval_poly(G(), newfe.ref_nodes[D], P)
    newns = nbr_ho_nodes(newfe)
    nv,nel = size(m.el)
    xdg = dg_nodes(m)
    newxdg = reshape(Pfe * (m.fe.shapefcn_coeff[D] * reshape(xdg, nv, D*nel)), newns*nel, D)
    neweldg = reshape(1:newns*nel, newns, nel)
    newx,newel = unique_mesh_nodes(newxdg, neweldg)
    HighOrderMesh{D,G,newp,T}(newfe, newx, newel, m.nbor)
end

change_ref_nodes(m::HighOrderMesh{D,G,P,T}, newsline::Vector{T}) where {D,G,P,T} =
    change_ref_nodes(m, FiniteElement(G(), newsline))

change_degree(m::HighOrderMesh{D,G,P,T}, newp::Int) where {D,G,P,T} =
    change_ref_nodes(m, FiniteElement(G(), newp, T))
