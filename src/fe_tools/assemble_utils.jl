using StaticArrays, LinearAlgebra, SparseArrays

struct FEM_precomp
    dim      # Dimension
    ng       # Nbr Gauss points
    ns       # Nbr solution points
    nel      # Nbr elements
    gξ       # Gauss nodes (ng x D)
    gw       # Gauss weights (ng)
    gϕ       # Basis functions at Gauss points (ng x ns)
    gϕξ      # Derivatives of gϕ w.r.t. ref coords (ng x ns x D)
    gJinv    # Element mapping inv Jacobian at all mesh Gauss points (ng x nel of SMatrix{D,D})
    gwJdet   # Element mapping Jacobian determinant at all mesh Gauss points (ng x nel)
    gx       # xyz-coordinates at all mesh Gauss points (ng x nel x D)
    gxξ      # Derivatives of gx w.r.t. ref coords (ng x nel x D x D)
end

function FEM_precomp(m::HighOrderMesh{D}; quadrature_degree=3*porder(m)) where {D}
    gξ,gw = quadrature(elgeom(m), quadrature_degree)
    ns,nel = size(m.el)
    ng = length(gw)

    gϕ = eval_shapefcns(m.fe, gξ, gradient=false)
    gϕξ = eval_shapefcns(m.fe, gξ, gradient=true)

    gx = eval_fcn(m.fe, dg_nodes(m), gξ, gradient=false)
    gxξ = eval_fcn(m.fe, dg_nodes(m), gξ, gradient=true)
    
    gJ = [ SMatrix{D,D}(view(gxξ,ig,iel,:,:)) for ig in 1:ng, iel in 1:nel ]
    gJinv = inv.(gJ)
    gwJdet = @. gw * det(gJ)

    FEM_precomp(D, ng, ns, nel, gξ, gw, gϕ, gϕξ, gJinv, gwJdet, gx, gxξ)
end

function eval_gϕx(pc::FEM_precomp, iel)
    gϕx = similar(pc.gϕξ)
    for ig = 1:pc.ng
        invJ = pc.gJinv[ig,iel]
        gϕx[ig,:,:] = pc.gϕξ[ig,:,:] * invJ
    end
    gϕx
end

function elmat_mass(pc::FEM_precomp, iel)
    cMel = pc.gϕ' * (pc.gwJdet[:,iel] .* pc.gϕ)
end

function elmat_laplace(pc::FEM_precomp, iel)
    gϕx = eval_gϕx(pc, iel)
    gwϕx = pc.gwJdet[:,iel] .* gϕx
    cAel = sum( (gϕx[:,:,d]' * gwϕx[:,:,d] for d = 1:pc.dim) )
end

function elres_rhs(pc::FEM_precomp, iel, fcn_rhs=x->1)
    f = [ fcn_rhs(view(pc.gx, ig, iel, :)) for ig = 1:pc.ng ]
    cfel = ((pc.gwJdet[:,iel] .* f)' * pc.gϕ)[:]
end

function assemble_matrix(el, fcn_Ael)
    ns,nel = size(el)
    ii = reshape(repeat(el, ns, 1), ns, ns, nel)
    jj = permutedims(ii, (2,1,3))
    Ael = cat( (fcn_Ael(iel) for iel = 1:nel)..., dims=3)
    sparse(ii[:], jj[:], Ael[:])
end

function assemble_vector(el, fcn_fel)
    ns,nel = size(el)
    fel = cat( (fcn_fel(iel) for iel = 1:nel)..., dims=2)
    f = zeros(eltype(fel), maximum(el))
    for iel in 1:nel
        f[el[:,iel]] .+= fel[:,iel]
    end
    f
end

function strong_dirichlet!(A, f, bndix)
    A[bndix,:] .= 0
    A[:,bndix] .= 0
    A[bndix,bndix] .= I(length(bndix))
    f[bndix] .= 0
end

function cg_mass(m::HighOrderMesh, pc::FEM_precomp)
    A = assemble_matrix(m.el, i -> elmat_mass(pc, i))
end


"""
    cg_poisson(m::HighOrderMesh, pc::FEM_precomp, fcn_rhs=xy->xy[1]^2, dirichlet_bnds=nothing)

Examples:
```julia
# Solve and plot -∇²u = x² with zero Dirichlet boundary conditions on the unit circle
m = mshcircle(2,p=3)
pc = FEM_precomp(m)
u,A,f = cg_poisson(m, pc, xy->xy[1]^2)
plot(m,u, contours=10, mesh_edges=true)
```
"""
function cg_poisson(m::HighOrderMesh, pc::FEM_precomp, fcn_rhs=xy->xy[1]^2, dirichlet_bnds=nothing)
    A = assemble_matrix(m.el, i -> elmat_laplace(pc, i))
    f = assemble_vector(m.el, i -> elres_rhs(pc, i, fcn_rhs))
    strong_dirichlet!(A, f, boundary_nodes(m, dirichlet_bnds))
    u = A \ f
    u,A,f
end

