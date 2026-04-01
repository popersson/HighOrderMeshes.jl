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
