###########################################################################
## DG operators: interior penalty Laplacian and upwind convection
#
# Both operators return the global sparse matrix on the element-local
# degrees of freedom, numbered i + ns(e-1) as in vec(U) for U of size
# ns × nel, and the right-hand side (ns × nel) contributed by the boundary
# data. Boundary conditions are a Dict from boundary numbers to functions
# g(x); boundaries not in the Dict get the natural condition: homogeneous
# Neumann for the Laplacian, outflow for convection.
#
# Face integrals are evaluated at the face quadrature points of the
# element with the lower number; the neighbor's traces are read at the
# matching points through nbpt. Face blocks are small BLAS products of
# trace matrices (ngf × ns) written into one preallocated block.

###########################################################################
## Block-sparse accumulation

# Storage for nblocks blocks of size ns × ns in coordinate format
mutable struct BlockCOO{T}
    ns::Int
    n::Int             # entries stored so far
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
end
BlockCOO(ns, nblocks, ::Type{T}=Float64) where {T} =
    BlockCOO{T}(ns, 0, Vector{Int}(undef, ns^2*nblocks), Vector{Int}(undef, ns^2*nblocks), Vector{T}(undef, ns^2*nblocks))

"""Add the block `A` (`ns × ns`) coupling test element `e1` with trial element `e2`."""
function addblock!(b::BlockCOO, A, e1, e2)
    ns, n = b.ns, b.n
    @inbounds for j in 1:ns, i in 1:ns
        n += 1
        b.I[n] = i + ns*(e1-1)
        b.J[n] = j + ns*(e2-1)
        b.V[n] = A[i, j]
    end
    b.n = n
end

SparseArrays.sparse(b::BlockCOO, ne) =
    sparse(view(b.I, 1:b.n), view(b.J, 1:b.n), view(b.V, 1:b.n), b.ns*ne, b.ns*ne)

"""Number of blocks of a DG matrix: one per element, four per interior face, one per boundary face satisfying `bndsel`."""
function nblocks(nb::AbstractMatrix{Neighbor}, bndsel)
    nel = size(nb, 2)
    nel + 2*count(n -> !isboundary(n), nb) + count(n -> isboundary(n) && bndsel(n), nb)
end

###########################################################################
## Traces on a face
##
## pts is the order of the face points: 1:ngf for the element whose face is
## integrated, nbpt[:,j,e] for its neighbor. wn is w·dS·n of the integrated
## side, so that Dn is the normal derivative in the direction of that n.

"""`T[g,i] = ϕ_i` of element `e` at the points `pts` of its local face `j`."""
function trace!(T, d::DGData, e, j, pts)
    r = d.fref[j]
    @inbounds for i in axes(T, 2), (g, k) in enumerate(pts)
        T[g, i] = r.ϕ[k, i]
    end
end

"""`Dn[g,i] = w·dS ∂ϕ_i/∂n` of element `e` at the points `pts` of its local face `j`."""
function ntrace!(Dn, d::DGData{D}, e, j, pts, wn) where {D}
    r   = d.fref[j]
    ngf = npoints(r)
    @inbounds for (g, k) in enumerate(pts)
        c = d.fmet.Jinv[k, j, e] * wn[g]                   # ∂ϕ/∂n w·dS = ∇ξϕ · (J⁻¹ w·dS·n)
        for i in axes(Dn, 2)
            Dn[g, i] = dot(c, gradvec(r.ϕξ, k, i, ngf, Val(D)))
        end
    end
end

###########################################################################
## Interior penalty Laplacian

"""
    dg_laplace(d::DGData, dirichlet) -> (A, b)

Symmetric interior penalty discretization of `-∇²u`: the matrix `A` and the
right-hand side `b` contributed by the Dirichlet data. `dirichlet` maps
boundary numbers to functions `g(x)`; other boundaries get a homogeneous
Neumann condition.

The penalty is the explicit expression of Shahbazi (2005): with `V_K` the
volume of element `K` and `S_K` its surface, counting Dirichlet faces twice
and Neumann faces not at all, `σ = c_p (S_K⁺/V_K⁺ + S_K⁻/V_K⁻)/4` on interior
faces and `σ = c_p S_K/(2V_K)` on Dirichlet faces, with the trace constant
`c_p = (p+1)(p+D)/D` on simplices and `(p+1)²` on blocks.
"""
function dg_laplace(d::DGData{D}, dirichlet) where {D}
    m, fm = d.m, d.fmet
    ns, ne, nf = nnodes(m.fe), nel(m), size(m.nb, 1)
    ngf = npoints(d.fref[1])
    T   = eltype(d.ref.ϕ)
    p   = porder(m)
    isdirichlet(nb) = isboundary(nb) && haskey(dirichlet, bndtag(nb))

    # penalty: c_p S_K / V_K per element
    cp   = elgeom(m) isa Simplex ? (p+1)*(p+D)/D : (p+1)^2
    vol  = vec(sum(d.met.wJ, dims=1))
    area = [ sum(norm, view(fm.wn, :, j, e)) for j in 1:nf, e in 1:ne ]
    surf = [ sum(area[j, e] * (isdirichlet(m.nb[j, e]) ? 2 : isboundary(m.nb[j, e]) ? 0 : 1) for j in 1:nf) for e in 1:ne ]
    σK   = cp .* surf ./ vol

    blocks = BlockCOO(ns, nblocks(m.nb, isdirichlet), T)
    b      = zeros(T, ns, ne)

    # volume terms
    K = elmats(elmat_laplace!, d)
    for e in 1:ne
        addblock!(blocks, view(K, :, :, e), e, e)
    end

    # face terms
    Tp, Tn, Dp, Dn, WT = (zeros(T, ngf, ns) for _ in 1:5)
    wdS, gq = zeros(T, ngf), zeros(T, ngf)
    blk = zeros(T, ns, ns)
    for e in 1:ne, j in 1:nf
        nb = m.nb[j, e]
        wn = view(fm.wn, :, j, e)
        wdS .= norm.(wn)
        if !isboundary(nb)
            nb.el > e || continue                                 # every interior face once
            trace!(Tp, d, e, j, 1:ngf);  ntrace!(Dp, d, e, j, 1:ngf, wn)
            pts = view(d.nbpt, :, j, e)
            trace!(Tn, d, nb.el, nb.face, pts);  ntrace!(Dn, d, nb.el, nb.face, pts, wn)
            σ = (σK[e] + σK[nb.el]) / 4
            # [u] = u⁺ - u⁻ with n pointing from + (element e) to - (the neighbor)
            sides = ((Tp, Dp, 1, e), (Tn, Dn, -1, Int(nb.el)))
            for (Tα, Dα, sα, eα) in sides, (Tβ, Dβ, sβ, eβ) in sides
                WT .= wdS .* Tβ
                mul!(blk, Tα', WT, σ*sα*sβ, false)                # σ [u][v]
                mul!(blk, Tα', Dβ, -sα/2, true)                   # -{∂ₙu}[v]
                mul!(blk, Dα', Tβ, -sβ/2, true)                   # -{∂ₙv}[u]
                addblock!(blocks, blk, eα, eβ)
            end
        elseif isdirichlet(nb)
            trace!(Tp, d, e, j, 1:ngf);  ntrace!(Dp, d, e, j, 1:ngf, wn)
            σ = σK[e] / 2
            WT .= wdS .* Tp
            mul!(blk, Tp', WT, σ, false)
            mul!(blk, Tp', Dp, -1, true)
            mul!(blk, Dp', Tp, -1, true)
            addblock!(blocks, blk, e, e)
            gq .= dirichlet[bndtag(nb)].(view(fm.x, :, j, e))     # σ ∫ g v - ∫ g ∂ₙv
            mul!(view(b, :, e), WT', gq, σ, true)
            mul!(view(b, :, e), Dp', gq, -1, true)
        end
    end
    sparse(blocks, ne), b
end

###########################################################################
## Upwind convection

"""
    dg_convection(d::DGData, vel, dirichlet) -> (A, b)

Upwind DG discretization of `∇·(v u)` with the velocity field `vel(x)`
(returning a D-vector): the matrix `A` and the right-hand side `b`
contributed by the inflow data. `dirichlet` maps boundary numbers to
functions `g(x)` giving `u` where the flow enters; outflow parts and
boundaries not in the Dict use the interior value.
"""
function dg_convection(d::DGData{D}, vel, dirichlet) where {D}
    m, r, met, fm = d.m, d.ref, d.met, d.fmet
    ns, ne, nf = nnodes(m.fe), nel(m), size(m.nb, 1)
    ng, ngf = npoints(r), npoints(d.fref[1])
    T = eltype(r.ϕ)
    blocks = BlockCOO(ns, nblocks(m.nb, nb -> true), T)
    b      = zeros(T, ns, ne)
    blk    = zeros(T, ns, ns)

    # volume terms: -∫ u v·∇ψ = -E'ϕ with E[g,i] = wJ (J⁻¹v)·∇ξψ_i
    E = zeros(T, ng, ns)
    for e in 1:ne
        @inbounds for g in 1:ng
            c = met.wJ[g, e] * (met.Jinv[g, e] * SVector{D}(vel(met.x[g, e])))
            for i in 1:ns
                E[g, i] = dot(c, gradvec(r.ϕξ, g, i, ng, Val(D)))
            end
        end
        mul!(blk, E', r.ϕ, -1, false)
        addblock!(blocks, blk, e, e)
    end

    # face terms: upwind flux w·dS (v·n)⁺ u⁺ + w·dS (v·n)⁻ u⁻ (positive and negative parts)
    Tp, Tn, F = (zeros(T, ngf, ns) for _ in 1:3)
    vn, gq = zeros(T, ngf), zeros(T, ngf)
    for e in 1:ne, j in 1:nf
        nb = m.nb[j, e]
        @inbounds for g in 1:ngf
            vn[g] = dot(SVector{D}(vel(fm.x[g, j, e])), fm.wn[g, j, e])
        end
        trace!(Tp, d, e, j, 1:ngf)
        if !isboundary(nb)
            nb.el > e || continue                                 # every interior face once
            trace!(Tn, d, nb.el, nb.face, view(d.nbpt, :, j, e))
            sides = ((Tp, 1, e), (Tn, -1, Int(nb.el)))            # (trace, sign of n, element)
            for (Tα, sα, eα) in sides, (Tβ, sβ, eβ) in sides
                F .= (sβ > 0 ? max.(vn, 0) : min.(vn, 0)) .* Tβ   # upwind: u⁺ where v·n > 0, u⁻ where v·n < 0
                mul!(blk, Tα', F, sα, false)
                addblock!(blocks, blk, eα, eβ)
            end
        else
            g = get(dirichlet, bndtag(nb), nothing)
            if g === nothing
                F .= vn .* Tp                                     # outflow: interior value
            else
                F .= max.(vn, 0) .* Tp                            # interior value where the flow leaves,
                gq .= g.(view(fm.x, :, j, e))                     # the data where it enters
                mul!(view(b, :, e), Tp', min.(vn, 0) .* gq, -1, true)
            end
            mul!(blk, Tp', F, true, false)
            addblock!(blocks, blk, e, e)
        end
    end
    sparse(blocks, ne), b
end
