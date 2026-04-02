###########################################################################
## Face-to-node mapping

"""
    mkface2nodes(eg::ElementGeometry, sface, svol)
    mkface2nodes(fe::FiniteElement)
    mkface2nodes(m::HighOrderMesh)

Build the face-to-node index map for element geometry `eg`. Returns an
`(nfacenodes × nfaces)` integer matrix where column `j` contains the local
node indices (within the element) that lie on face `j`.

Matching is done by comparing shape function values at face reference nodes
against volume reference nodes, so it is basis-order-aware.
"""
function mkface2nodes(eg::ElementGeometry{D}, sface::AbstractArray{T}, svol::AbstractArray{T}) where {D,T}
    fmap       = facemap(eg)
    basis_face = snap.(eval_shapefcns(subgeom(eg, D-1), sface))
    basis_vol  = snap.(eval_shapefcns(eg, svol))
    f2n = fill(0, size(basis_face,1), size(fmap,2))
    for (ii, ic) in enumerate(eachcol(fmap))
        f2n[:,ii] .= indexin(eachrow(basis_face), eachrow(basis_vol[:,ic]))
    end
    f2n
end

mkface2nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T} =
    mkface2nodes(G(), ref_nodes(fe, D-1), ref_nodes(fe, D))

mkface2nodes(m::HighOrderMesh) = mkface2nodes(m.fe)

###########################################################################
## LDG switch functions

"""
    mkldgswitch(eg::Block, nb) -> sw
    mkldgswitch(m::HighOrderMesh) -> sw

Compute the LDG (Local Discontinuous Galerkin) switch matrix for a block mesh.
Returns an `(nfaces × nelems)` integer matrix `sw` with values in `{0, 1}`,
where `sw[j, iel]` is the switch function on face `j` of element `iel`.
Adjacent elements always have opposite switch values on their shared face,
ensuring a consistent upwind direction for the auxiliary variable flux.

Only implemented for `Block` geometries.
"""
function mkldgswitch(eg::Block{D}, nb) where {D}
    nf, nel = size(nb)
    sw = fill(-1, nf, nel)  # -1 = not yet assigned

    # On a structured block mesh, faces come in opposite pairs (e.g. left/right).
    opposite_face(i) = i - iseven(i) + isodd(i)

    while true
        idx = findfirst(isequal(-1), sw)
        isnothing(idx) && break
        j, iel = Tuple(idx)

        # Walk along a line of elements in one direction, then back the other way
        iel0, j0, dir = iel, j, 1
        while true
            sw[j,iel] = dir
            jel, k, _ = nb[j,iel]
            if jel > 0
                sw[k,jel] = 1 - dir  # neighbor gets opposite switch
                iel = jel
                j   = opposite_face(k)
                if sw[j,iel] != -1
                    if dir == 1
                        dir = 0;  iel = iel0;  j = opposite_face(j0)  # reverse
                    else
                        break
                    end
                end
            elseif dir == 1
                dir = 0;  iel = iel0;  j = opposite_face(j0)  # hit boundary, reverse
            else
                break
            end
        end
    end
    return sw
end

mkldgswitch(m::HighOrderMesh) = mkldgswitch(elgeom(m), m.nb)

###########################################################################
## Node/face reordering for LDG compatibility

"""
    align_with_ldgswitch!(m::HighOrderMesh{2,Block{2},P}, sw=nothing)

Reorder element nodes and neighbor data in-place so that the local node
numbering of each 2D quad element is consistent with the LDG switch `sw`.
Computes `sw` via `mkldgswitch` if not supplied.

The four possible switch patterns `[s1,s2,s3,s4] ∈ {0,1}^4` on the four faces
correspond to four 90° rotations of the reference element. This function
applies the appropriate rotation to `m.el` and updates `m.nb` accordingly.
"""
function align_with_ldgswitch!(m::HighOrderMesh{2,Block{2},P}, sw=nothing) where P
    sw = isnothing(sw) ? mkldgswitch(m) : sw

    # The four distinct switch patterns on a quad's faces, one per rotation
    switch_cases = [ [1,0,1,0], [1,0,0,1], [0,1,0,1], [0,1,1,0] ]
    mapcase = [ findfirst([csw] .== switch_cases) for csw in eachcol(sw) ]

    # Node reordering maps: successive 90° rotations of the (P+1)×(P+1) grid
    nds    = collect(reshape(1:(P+1)^2, P+1, P+1))
    ndmaps = [nds]
    for _ = 2:4
        push!(ndmaps, rotl90(ndmaps[end]))
    end
    ndmaps = reshape.(ndmaps, :)

    # Face reordering maps and their inverses (one per rotation)
    fcmaps  = [ [1,2,3,4], [4,3,1,2], [2,1,4,3], [3,4,2,1] ]
    ifcmaps = [ invperm(fc) for fc in fcmaps ]

    for iel = 1:size(m.el,2)
        cmap = mapcase[iel]
        m.el[:,iel] .= m.el[:,iel][ndmaps[cmap]]
        cnb = m.nb[:,iel]
        for j = 1:4
            if cnb[j][2] > 0
                jel    = cnb[j][1]
                cnb[j] = (jel, ifcmaps[mapcase[jel]][cnb[j][2]], 0)
            end
        end
        m.nb[:,iel] = cnb[fcmaps[cmap]]
    end
end
