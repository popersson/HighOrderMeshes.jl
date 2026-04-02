###########################################################################
## Legendre polynomials

"""
    legendre_poly(x, p; diff=false)

Evaluate Legendre polynomials `L_0, ..., L_p` at each point in `x` on `[-1,1]`
using the three-term recurrence. Returns an `(length(x) × p+1)` matrix.
If `diff=true`, returns `(y, dy)` where `dy` holds the derivatives.
"""
function legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T
    y = ones(T, length(x), p+1)
    y[:,2] .= x
    if diff
        dy = zeros(T, length(x), p+1)
        dy[:,2] .= one(T)
    end
    for j = 1:p-1
        @. @views y[:,j+2]  = ((2j+1)*x*y[:,j+1]  - j*y[:,j])  / (j+1)
        if diff
            @. @views dy[:,j+2] = ((2j+1)*(x*dy[:,j+1] + y[:,j+1]) - j*dy[:,j]) / (j+1)
        end
    end
    diff ? (y, dy) : y
end

"""
    legendre01_poly(x, p; diff=false)

Shifted Legendre polynomials on `[0,1]`, obtained by the substitution `x → 2x-1`.
Same return convention as `legendre_poly`.
"""
function legendre01_poly(x::AbstractVector{T}, p::Int; diff=false) where T
    ys = legendre_poly(2x .- 1, p; diff=diff)
    diff && (ys[2] .*= 2)  # chain rule: d/dx = 2 * d/d(2x-1)
    ys
end

###########################################################################
## Multivariate polynomial bases

"""
    multivar_legendre01_poly(x, p; gradient=false)

Tensor-product shifted Legendre basis of degree `p` on `[0,1]^D`,
where `D = size(x,2)` and `x` is `npts × D`.

- `gradient=false` → `npts × (p+1)^D`
- `gradient=true`  → `npts × (p+1)^D × D`

Used as the polynomial basis for `Block` elements.
"""
function multivar_legendre01_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
    D = size(x,2)
    ix = [ i for i in Iterators.product(fill(1:p+1,D)...)][:]
    # Evaluate 1D Legendre polynomials along each coordinate axis
    Ls = [ legendre01_poly(xx, p) for xx in eachcol(x) ]

    # Build tensor-product values by multiplying 1D factors
    function make_outer_prod!(Ls, P)
        for i = 1:size(x,1), j = 1:length(ix)
            val = one(T)
            for k = 1:D
                val *= Ls[k][i, ix[j][k]]
            end
            P[i,j] = val
        end
    end

    if !gradient
        P = zeros(T, size(x,1), length(ix))
        make_outer_prod!(Ls, P)
    else
        # Gradient: replace the k-th factor with its derivative (product rule)
        dLs = [ legendre01_poly(xx, p, diff=true)[2] for xx in eachcol(x) ]
        P = zeros(T, size(x,1), length(ix), D)
        for k = 1:D
            cLs = [ i == k ? dLs[i] : Ls[i] for i = 1:D ]
            make_outer_prod!(cLs, view(P, :, :, k))
        end
    end
    P
end

"""
    multivar_monomial_poly(x, p; gradient=false)

Monomial basis `{x₁^i₁ ⋯ xD^iD : i₁+⋯+iD ≤ p}` on `ℝ^D`,
where `D = size(x,2)` and `x` is `npts × D`. Monomials are ordered by
total degree.

- `gradient=false` → `npts × nbasis`
- `gradient=true`  → `npts × nbasis × D`

Used as the polynomial basis for `Simplex` elements.
"""
function multivar_monomial_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
    isempty(x) && return [one(T);;]
    x = x'
    D = size(x,1)
    powers = sort([ i for i in Iterators.product(fill(0:p,D)...) if sum(i) <= p ],
                  lt=(i,j) -> sum(i) < sum(j))
    if !gradient
        return [ prod(xx .^ k) for xx in eachcol(x), k in powers ]
    else
        # Differentiate each monomial: d/dx_d (x^k) = k_d * x^(k - e_d)
        q = zeros(T, size(x,2), length(powers), D)
        for d = 1:D
            mul        = (ix -> ix[d]).(powers)
            mask       = Tuple(Int(i==d) for i = 1:D)
            diffpowers = [ max.(0, power .- mask) for power in powers ]
            q[:,:,d]   = [ m .* prod(xx .^ k) for xx in eachcol(x), (k,m) in zip(diffpowers,mul) ]
        end
        q
    end
end

"""
    eval_poly(eg::ElementGeometry, s, p; gradient=false)

Evaluate the polynomial basis for geometry `eg` at reference points `s` up to
degree `p`. Dispatches to `multivar_monomial_poly` for `Simplex` and
`multivar_legendre01_poly` for `Block`.
"""
eval_poly(eg::Simplex, s, p; gradient=false) = multivar_monomial_poly(s, p; gradient=gradient)
eval_poly(eg::Block,   s, p; gradient=false) = multivar_legendre01_poly(s, p; gradient=gradient)

###########################################################################
## 1D quadrature nodes and weights
#
# Nodes are computed by Newton iteration on the Legendre polynomial roots.
# The [0,1] variants are simple affine rescalings of the [-1,1] versions.

"""
    gauss_legendre_nodes(n; T=Float64)
    gauss_legendre01_nodes(n; T=Float64)

Return `n` Gauss-Legendre nodes on `[-1,1]` (or `[0,1]`), computed by Newton
iteration on `L_n`. These are the optimal interior quadrature nodes for
polynomials of degree up to `2n-1`.
"""
function gauss_legendre_nodes(n::Int; T::Type=Float64)
    x = [ cos(T(π) / 4n * (4k - 2)) for k = n:-1:1 ]  # initial guess
    for _ = 1:100
        P, dP = legendre_poly(x, n, diff=true)
        xold  = copy(x)
        @. x  = xold - P[:,end] / dP[:,end]
        maximum(abs.(x - xold)) < 2*eps(T) && return x
    end
    error("No convergence in Gauss-Legendre Newton iterations")
end

gauss_legendre01_nodes(n::Int; T::Type=Float64) = (gauss_legendre_nodes(n; T) .+ 1) ./ 2

"""
    gauss_legendre_quadrature(n; T=Float64)
    gauss_legendre01_quadrature(n; T=Float64)

Gauss-Legendre nodes and weights on `[-1,1]` (or `[0,1]`) with `n` points.
Degree of precision: `2n-1`.

```julia
x, w = gauss_legendre_quadrature(3)
w' * x.^4   # ≈ 2/5  (exact: DoP=5 ≥ 4)
```
"""
function gauss_legendre_quadrature(n::Int; T::Type=Float64)
    x    = gauss_legendre_nodes(n; T)
    P,dP = legendre_poly(x, n, diff=true)
    w    = @. 2 / (1 - x^2) / dP[:,end]^2
    x, w
end

function gauss_legendre01_quadrature(n::Int; T::Type=Float64)
    x, w = gauss_legendre_quadrature(n; T)
    (x .+ 1) ./ 2, w ./ 2
end

"""
    gauss_lobatto_nodes(n; T=Float64)
    gauss_lobatto01_nodes(n; T=Float64)

Return `n` Gauss-Lobatto nodes on `[-1,1]` (or `[0,1]`). Includes the
endpoints `±1`; interior nodes are computed by Newton iteration.
"""
function gauss_lobatto_nodes(n::Int; T::Type=Float64)
    n1 = n - 1
    x  = @. cos(π * ((n1:-1:0) / T(n1)))  # initial guess (Chebyshev)
    for _ = 1:100
        P    = legendre_poly(x, n1)
        xold = copy(x)
        x    = @. xold - (x * P[:,end] - P[:,end-1]) / ((n1+1) * P[:,end])
        maximum(abs.(x - xold)) < 2*eps(T) && return x
    end
    error("No convergence in Gauss-Lobatto Newton iterations")
end

gauss_lobatto01_nodes(n::Int; T::Type=Float64) = (gauss_lobatto_nodes(n; T) .+ 1) ./ 2

"""
    gauss_lobatto_quadrature(n; T=Float64)
    gauss_lobatto01_quadrature(n; T=Float64)

Gauss-Lobatto nodes and weights on `[-1,1]` (or `[0,1]`) with `n` points.
Degree of precision: `2n-3`.

```julia
x, w = gauss_lobatto_quadrature(4)
w' * x.^4   # ≈ 2/5  (exact: DoP=5 ≥ 4)
```
"""
function gauss_lobatto_quadrature(n::Int; T::Type=Float64)
    x = gauss_lobatto_nodes(n; T)
    P = legendre_poly(x, n-1)
    w = @. 2 / ((n-1) * n * P[:,end]^2)
    x, w
end

function gauss_lobatto01_quadrature(n::Int; T::Type=Float64)
    x, w = gauss_lobatto_quadrature(n; T)
    (x .+ 1) ./ 2, w ./ 2
end

###########################################################################
## Reference node distributions

"""Equispaced nodes on `[0,1]` at polynomial order `p` (as exact rationals)."""
equispaced(p::Int) = (0:p) // p

"""
    ref_nodes(::Block{D}, s1)
    ref_nodes(::Simplex{D}, s1)

Build `D`-dimensional reference nodes from a 1D node vector `s1` on `[0,1]`.

- `Block`: full tensor product `s1^D`, giving `(p+1)^D` nodes.
- `Simplex`: tensor product filtered to `sum(x) ≤ 1`, giving `binomial(p+D,D)` nodes.
"""
function ref_nodes(::Block{D}, s1::AbstractVector{T}) where {D,T}
    sD = zeros(T, 1, 0)
    for d = 1:D
        sD = hcat(repeat(sD, length(s1), 1),
                  reshape(repeat(s1', size(sD,1), 1), :, 1))
    end
    sD
end

function ref_nodes(::Simplex{D}, s1) where {D}
    sD = ref_nodes(Block{D}(), s1)
    sD[sum(sD, dims=2)[:,1] .<= 1, :]
end

# Identify element geometry from dimension D, polynomial order P, and node count.
# Companion to find_elgeom(D, nv) in element_geometry.jl, which infers from vertex count.
function find_elgeom(D, P, nnodes)
    (P+1)^D       == nnodes && return Block{D}()
    binomial(P+D,D) == nnodes && return Simplex{D}()
    error("Cannot determine element geometry for D=$D, P=$P, nnodes=$nnodes")
end

###########################################################################
## Multi-dimensional quadrature

"""
    quadrature(eg::ElementGeometry, p; T=Float64)

Return `(ξ, w)`: quadrature nodes `ξ` (`npts × D`) and weights `w` (`npts`)
exact for polynomials of degree `p` on element geometry `eg`.

- `Block`: tensor-product Gauss-Legendre rule on `[0,1]^D`.
- `Simplex`: hard-coded symmetric rules from `simplex_quadrature.jl`.
"""
quadrature(eg::ElementGeometry, p) = error("Quadrature not implemented for $(typeof(eg))")

quadrature(eg::Simplex{D}, p) where {D} = simplex_quadrature(SimplexQuadRule{D,p}())

function quadrature(eg::Block{D}, p; T::Type=Float64) where D
    ξ0, w0 = gauss_legendre01_quadrature(p; T)
    ξ = ref_nodes(eg, ξ0)
    # Tensor-product weights: outer product of 1D weight vectors
    w = w0
    for d = 2:D
        w = (w .* w0')[:]
    end
    ξ, w
end

