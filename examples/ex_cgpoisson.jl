###
### Poisson in 2D using CG finite elements
###

using HighOrderMeshes
using Plots

using StaticArrays, LinearAlgebra, SparseArrays

function FEMPrecomp(m::HighOrderMesh; quadrature_degree=3*porder(m))
    gx,gw = quadrature(elgeom(m), quadrature_degree)
    ϕ = eval_shapefcns(m.fe, gx, gradient=false)
    ϕξ = eval_shapefcns(m.fe, gx, gradient=true)
    (gx=gx, gw, ϕ, ϕξ)
end

function cg_poisson(m::HighOrderMesh, pc)
    gJ = eval_fcn(m.fe, dg_nodes(m), pc.gx, gradient=true)

    ns = nbr_ho_nodes(m.fe)
    nel = size(m.el,2)
    ng = size(pc.gx,1)

    Ael = zeros(ns,ns,nel)
    fel = zeros(ns,nel)

    for iel in 1:nel, ig in 1:ng
        J = SMatrix{2,2}(view(gJ,ig,iel,:,:))
        detJ = det(J)
        invJ = inv(J)
        
        ϕx = pc.ϕξ[ig,:,:] * invJ
        Ael[:,:,iel] .+= pc.gw[ig] * detJ * (ϕx[:,1]*ϕx[:,1]' + ϕx[:,2]*ϕx[:,2]')
        fel[:,iel] .+= pc.gw[ig] * detJ * pc.ϕ[ig,:]
    end
    Ael,fel
end

function assemble_matrix(el, Ael)
    ns,nel = size(el)
    ii = reshape(repeat(el, ns, 1), ns, ns, nel)
    jj = permutedims(ii, (2,1,3))
    sparse(ii[:], jj[:], Ael[:])
end

function assemble_vector(el, fel)
    ns,nel = size(el)
    f = zeros(eltype(fel), maximum(el))
    for iel in 1:nel
        f[el[:,iel]] .+= fel[:,iel]
    end
    f
end

function strong_dirichlet(A, f, bndix)
    A[bndix,:] .= 0
    A[:,bndix] .= 0
    A[bndix,bndix] .= I(length(bndix))
    f[bndix] .= 0
    A,f
end



m = mshcircle(2,p=3)
pc = FEMPrecomp(m)
Ael,fel = cg_poisson(m, pc)
A = assemble_matrix(m.el, Ael)
f = assemble_vector(m.el, fel)

A,f = strong_dirichlet(A, f, boundary_nodes(m,1:4))
u = A \ f
plot(m,u, contours=10, mesh_edges=true)

