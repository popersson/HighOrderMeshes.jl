###########################################################################
## Elemental kernels and continuous Galerkin assembly
#
# Elemental matrices are written into slice e of a preallocated ns × ns × nel
# array; B is a workspace of (at least) ng·D × ns. All matrix products are
# BLAS calls on views; the per-point loops apply the metric terms.

"""
    elmat_laplace!(K, d::FEMData, e, B)

Elemental Laplace matrix `K[:,:,e][i,j] = ∫ ∇ϕ_i·∇ϕ_j` of element `e`, computed
as `ϕξ' B` with `B = (wJ J⁻¹J⁻ᵀ) ∇ξϕ` at every quadrature point. `B` is a
workspace of size `ng·D × ns`.
"""
function elmat_laplace!(K, d::FEMData{D}, e, B) where {D}
    r, met = d.ref, d.met
    ng, ns = npoints(r), size(r.ϕ, 2)
    @inbounds for g in 1:ng
        Ji = met.Jinv[g, e]
        G  = met.wJ[g, e] * (Ji * Ji')
        for i in 1:ns
            setgradvec!(B, G * gradvec(r.ϕξ, g, i, ng, Val(D)), g, i, ng)
        end
    end
    mul!(view(K, :, :, e), r.ϕξ', view(B, 1:ng*D, :))
end

"""
    elmat_mass!(M, d::FEMData, e, B)

Elemental mass matrix `M[:,:,e] = ϕ' (wJ .* ϕ)` of element `e`. `B` is a
workspace with at least `ng` rows and `ns` columns.
"""
function elmat_mass!(M, d::FEMData, e, B)
    r, met = d.ref, d.met
    ng = npoints(r)
    Bm = view(B, 1:ng, :)
    Bm .= view(met.wJ, :, e) .* r.ϕ
    mul!(view(M, :, :, e), r.ϕ', Bm)
end

"""
    elmats(kernel!, d::FEMData) -> K

All elemental matrices of `kernel!` (for example [`elmat_laplace!`](@ref) or
[`elmat_mass!`](@ref)) as an `ns × ns × nel` array.
"""
function elmats(kernel!, d::FEMData{D}) where {D}
    ng, ns, ne = npoints(d.ref), nnodes(d.m.fe), nel(d.m)
    K = zeros(eltype(d.ref.ϕ), ns, ns, ne)
    B = zeros(eltype(d.ref.ϕ), ng*D, ns)
    for e in 1:ne
        kernel!(K, d, e, B)
    end
    K
end

"""
    elvec_source(d::FEMData, f) -> F

Elemental source vectors `F[i,e] = ∫ f ϕ_i` over element `e` for the function
`f(x)`, as an `ns × nel` array.
"""
elvec_source(d::FEMData, f) = d.ref.ϕ' * (d.met.wJ .* f.(d.met.x))

"""
    laplace_residual!(R, Q, d::FEMData, U)

Matrix-free `R = K U` for all elements (`U`, `R`: `ns × nel`) without forming
the elemental Laplace matrices: `Q = ∇ξϕ U`, then `Q ← wJ J⁻¹J⁻ᵀ Q` at every
point, then `R = ∇ξϕ' Q`. `Q` is a workspace of size `ng·D × nel`.
"""
function laplace_residual!(R, Q, d::FEMData{D}, U) where {D}
    r, met = d.ref, d.met
    ng, ne = npoints(r), nel(d.m)
    mul!(Q, r.ϕξ, U)
    @inbounds for e in 1:ne, g in 1:ng
        Ji = met.Jinv[g, e]
        setgradvec!(Q, met.wJ[g, e] * (Ji * (Ji' * gradvec(Q, g, e, ng, Val(D)))), g, e, ng)
    end
    mul!(R, r.ϕξ', Q)
end

###########################################################################
## Continuous Galerkin assembly

"""
    assemble_matrix(el, K) -> A

Global sparse matrix on the mesh nodes from the elemental matrices `K`
(`ns × ns × nel`) and the connectivity `el`.
"""
function assemble_matrix(el, K::AbstractArray{<:Any,3})
    ns, ne = size(el)
    ii = reshape(repeat(el, ns, 1), ns, ns, ne)      # ii[i,j,e] = el[i,e]
    jj = permutedims(ii, (2, 1, 3))                  # jj[i,j,e] = el[j,e]
    sparse(vec(ii), vec(jj), vec(K), maximum(el), maximum(el))
end

"""
    assemble_vector(el, F) -> f

Global vector on the mesh nodes from the elemental vectors `F` (`ns × nel`)
and the connectivity `el`.
"""
function assemble_vector(el, F::AbstractMatrix)
    f = zeros(eltype(F), maximum(el))
    for e in axes(el, 2)
        f[view(el, :, e)] .+= view(F, :, e)
    end
    f
end

"""
    strong_dirichlet(A, f, nodes, ud=0) -> (A, f)

Impose `u = ud` strongly at the given nodes of the system `A u = f`: the
rows and columns of `A` are replaced by the identity and the known values are
moved to the right-hand side, which keeps a symmetric `A` symmetric. `ud` is
a vector of values at `nodes`, or zero.
"""
function strong_dirichlet(A, f, nodes, ud=zeros(eltype(f), length(nodes)))
    n    = length(f)
    keep = ones(eltype(f), n); keep[nodes] .= 0
    u    = zeros(eltype(f), n); u[nodes] = ud
    P, Pd = Diagonal(keep), Diagonal(1 .- keep)
    P * A * P + Pd, P * (f - A * u) + u
end
