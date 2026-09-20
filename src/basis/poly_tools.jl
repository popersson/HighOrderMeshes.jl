###########################################################################
## Jacobi and Legendre polynomials
#
# All polynomial evaluations are scalar functions of a scalar x. Evaluate at
# many points by broadcasting or comprehensions. Exact number types such as
# Rational are preserved; integers are promoted through division.

"""
    jacobi(n, α, β, x)

Jacobi polynomial ``P_n^{(α,β)}(x)`` at the scalar `x`, evaluated by the
three-term recurrence. Vectors are handled by broadcasting, `jacobi.(n, α, β, x)`.
"""
function jacobi(n::Integer, α::Real, β::Real, x::Real)
    x  = x / one(x)                       # promote integers to a division-closed type
    n  < 0 && return zero(x)
    n == 0 && return one(x)
    p0 = one(x)
    p1 = ((α + β + 2) * x + (α - β)) / 2
    for i in 1:n-1
        a1 = 2 * (i + 1) * (i + α + β + 1) * (2i + α + β)
        a2 = (2i + α + β + 1) * (α^2 - β^2)
        a3 = (2i + α + β) * (2i + α + β + 1) * (2i + α + β + 2)
        a4 = 2 * (i + α) * (i + β) * (2i + α + β + 2)
        p0, p1 = p1, ((a2 + a3 * x) * p1 - a4 * p0) / a1
    end
    p1
end

"""
    djacobi(n, α, β, x)

Derivative of the Jacobi polynomial, using
``d/dx\\, P_n^{(α,β)}(x) = (n+α+β+1)/2 \\; P_{n-1}^{(α+1,β+1)}(x)``.
"""
djacobi(n::Integer, α::Real, β::Real, x::Real) =
    n <= 0 ? zero(x / one(x)) : (n + α + β + 1) * jacobi(n - 1, α + 1, β + 1, x) / 2

"""
    legendre(n, x),   dlegendre(n, x)
    legendre01(n, x), dlegendre01(n, x)

Legendre polynomial ``P_n(x)`` on `[-1,1]` and its derivative, and the shifted
versions on `[0,1]` obtained by the substitution `x → 2x-1`.

```julia
V = [legendre(n, x) for x in xs, n in 0:p]   # Vandermonde matrix
```
"""
legendre(n::Integer, x::Real)    = jacobi(n, 0, 0, x)
dlegendre(n::Integer, x::Real)   = djacobi(n, 0, 0, x)
legendre01(n::Integer, x::Real)  = legendre(n, 2x - 1)
dlegendre01(n::Integer, x::Real) = 2 * dlegendre(n, 2x - 1)

###########################################################################
## Multi-indices and node sets
#
# The canonical node order is the order of multiindices(eg, p): the first
# index varies fastest, and simplices keep only sum(i) <= p. The Gmsh, VTK and
# 3DG permutation tables are expressed relative to this order. Never change it.

"""
    multiindices(eg::ElementGeometry, p)

Multi-indices `(i_1, …, i_D)` of the degree-`p` polynomial space on `eg`, in
the canonical order: first index varying fastest, `0 ≤ i_d ≤ p`, and for
simplices `sum(i) ≤ p`.
"""
multiindices(::Block{D}, p::Integer) where {D} =
    vec(collect(Iterators.product(ntuple(_ -> 0:p, D)...)))
multiindices(::Simplex{D}, p::Integer) where {D} =
    filter(i -> sum(i; init=0) <= p, vec(collect(Iterators.product(ntuple(_ -> 0:p, D)...))))

"""Equispaced nodes on `[0,1]` at polynomial degree `p` (as exact rationals)."""
equispaced(p::Integer) = (0:p) // p

"""
    equispaced_nodes(eg::ElementGeometry, p)

Equispaced Lagrange nodes of degree `p` on the reference element, as an
`nnodes × D` matrix of rationals in the canonical node order.
"""
function equispaced_nodes(eg::ElementGeometry{D}, p::Integer) where {D}
    D == 0 && return zeros(Rational{Int}, 1, 0)
    p >= 1 || throw(ArgumentError("polynomial degree must be at least 1"))
    idx = multiindices(eg, p)
    [ i[d] // p for i in idx, d in 1:D ]
end

"""
    tensor_nodes(eg::Block{D}, s1)

Tensor product of the 1D node line `s1` on `[0,1]`, as a `length(s1)^D × D`
matrix in the canonical node order.
"""
function tensor_nodes(::Block{D}, s1::AbstractVector{T}) where {D,T}
    D == 0 && return zeros(T, 1, 0)
    idx = vec(collect(Iterators.product(ntuple(_ -> eachindex(s1), D)...)))
    [ s1[i[d]] for i in idx, d in 1:D ]
end

###########################################################################
## Polynomial bases on the reference elements

# Output number type for computations with points of element type T.
_outT(::AbstractArray{T}) where {T} = typeof(one(T) / one(T))

"""
    polybasis(eg::ElementGeometry, ξ, p)
    dpolybasis(eg::ElementGeometry, ξ, p)

Polynomial basis of degree `p` on the reference element `eg`, evaluated at the
points `ξ` (`nξ × D`). `polybasis` returns an `nξ × nbasis` matrix and
`dpolybasis` the gradients as an `nξ × nbasis × D` array. Basis functions are
ordered as `multiindices(eg, p)`.

- `Block`: tensor product of shifted Legendre polynomials on `[0,1]^D`.
- `Simplex`: orthogonal (unnormalized) Koornwinder basis in collapsed
  coordinates, following Hesthaven & Warburton.
"""
function polybasis(eg::Block{D}, ξ::AbstractVecOrMat, p::Integer) where {D}
    _checkdim(ξ, D)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    L   = [ T(legendre01(n, ξ[k, d])) for k in axes(ξ, 1), d in 1:D, n in 0:p ]
    [ prod(L[k, d, i[d]+1] for d in 1:D; init=one(T)) for k in axes(ξ, 1), i in idx ]
end

function dpolybasis(eg::Block{D}, ξ::AbstractVecOrMat, p::Integer) where {D}
    _checkdim(ξ, D)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    L   = [ T(legendre01(n, ξ[k, d]))  for k in axes(ξ, 1), d in 1:D, n in 0:p ]
    dL  = [ T(dlegendre01(n, ξ[k, d])) for k in axes(ξ, 1), d in 1:D, n in 0:p ]
    out = Array{T}(undef, size(ξ, 1), length(idx), D)
    for (j, i) in enumerate(idx), k in axes(ξ, 1), e in 1:D
        out[k, j, e] = prod(d == e ? dL[k, d, i[d]+1] : L[k, d, i[d]+1] for d in 1:D; init=one(T))
    end
    out
end

# Collapsed coordinates of the [0,1] triangle and tetrahedron. The singular
# vertex is mapped to a = -1 (b = -1); the value is irrelevant there because
# the affected basis functions carry a vanishing factor.
function _collapsed(x, y)
    a = y == 1 ? -one(x) : 2x / (1 - y) - 1
    b = 2y - 1
    a, b
end
function _collapsed(x, y, z)
    a = y + z == 1 ? -one(x) : 2x / (1 - y - z) - 1
    b = z == 1     ? -one(x) : 2y / (1 - z) - 1
    c = 2z - 1
    a, b, c
end

function polybasis(eg::Simplex{2}, ξ::AbstractVecOrMat, p::Integer)
    _checkdim(ξ, 2)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    V   = Matrix{T}(undef, size(ξ, 1), length(idx))
    for (n, (i, j)) in enumerate(idx), k in axes(ξ, 1)
        x, y = T(ξ[k, 1]), T(ξ[k, 2])
        a, b = _collapsed(x, y)
        V[k, n] = jacobi(i, 0, 0, a) * (1 - y)^i * jacobi(j, 2i + 1, 0, b)
    end
    V
end

function dpolybasis(eg::Simplex{2}, ξ::AbstractVecOrMat, p::Integer)
    _checkdim(ξ, 2)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    out = Array{T}(undef, size(ξ, 1), length(idx), 2)
    for (n, (i, j)) in enumerate(idx), k in axes(ξ, 1)
        x, y = T(ξ[k, 1]), T(ξ[k, 2])
        a, b = _collapsed(x, y)
        hb   = (1 - b) / 2                          # = 1 - y
        fa, dfa = jacobi(i, 0, 0, a),      djacobi(i, 0, 0, a)
        gb, dgb = jacobi(j, 2i + 1, 0, b), djacobi(j, 2i + 1, 0, b)
        hbm = i >= 1 ? hb^(i - 1) : one(T)         # (1-b)/2 to the power i-1
        dr  = dfa * gb * hbm
        tmp = dgb * hb^i
        i >= 1 && (tmp -= i * gb * hbm / 2)
        ds  = dfa * gb * (1 + a) / 2 * hbm + fa * tmp
        out[k, n, 1] = 2dr                          # d/dx = 2 d/dr on [0,1]
        out[k, n, 2] = 2ds
    end
    out
end

function polybasis(eg::Simplex{3}, ξ::AbstractVecOrMat, p::Integer)
    _checkdim(ξ, 3)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    V   = Matrix{T}(undef, size(ξ, 1), length(idx))
    for (n, (i, j, l)) in enumerate(idx), k in axes(ξ, 1)
        x, y, z = T(ξ[k, 1]), T(ξ[k, 2]), T(ξ[k, 3])
        a, b, c = _collapsed(x, y, z)
        hb = (1 - b) / 2
        V[k, n] = jacobi(i, 0, 0, a) * jacobi(j, 2i + 1, 0, b) * hb^i *
                  jacobi(l, 2i + 2j + 2, 0, c) * (1 - z)^(i + j)
    end
    V
end

function dpolybasis(eg::Simplex{3}, ξ::AbstractVecOrMat, p::Integer)
    _checkdim(ξ, 3)
    T   = _outT(ξ)
    idx = multiindices(eg, p)
    out = Array{T}(undef, size(ξ, 1), length(idx), 3)
    for (n, (i, j, l)) in enumerate(idx), k in axes(ξ, 1)
        x, y, z = T(ξ[k, 1]), T(ξ[k, 2]), T(ξ[k, 3])
        a, b, c = _collapsed(x, y, z)
        hb, hc  = (1 - b) / 2, (1 - c) / 2
        fa, dfa = jacobi(i, 0, 0, a),           djacobi(i, 0, 0, a)
        gb, dgb = jacobi(j, 2i + 1, 0, b),      djacobi(j, 2i + 1, 0, b)
        hl, dhl = jacobi(l, 2i + 2j + 2, 0, c), djacobi(l, 2i + 2j + 2, 0, c)
        hbm = i >= 1     ? hb^(i - 1)     : one(T)
        hcm = i + j >= 1 ? hc^(i + j - 1) : one(T)
        # r-derivative
        vr = dfa * gb * hl * hbm * hcm
        # s-derivative
        tmp = dgb * hb^i
        i >= 1 && (tmp -= i * gb * hbm / 2)
        tmp = fa * tmp * hl * hcm
        vs  = (1 + a) / 2 * vr + tmp
        # t-derivative
        vt  = (1 + a) / 2 * vr + (1 + b) / 2 * tmp
        tmp = dhl * hc^(i + j)
        i + j >= 1 && (tmp -= (i + j) * hl * hcm / 2)
        vt += fa * gb * tmp * hb^i
        out[k, n, 1] = 2vr
        out[k, n, 2] = 2vs
        out[k, n, 3] = 2vt
    end
    out
end

_checkdim(ξ, D) = size(ξ, 2) == D ||
    throw(DimensionMismatch("expected $D reference coordinates per point, got $(size(ξ, 2))"))

###########################################################################
## 1D quadrature nodes and weights
#
# Nodes are computed by Newton iteration on Legendre polynomial roots. The
# [0,1] variants are affine rescalings of the [-1,1] versions.

"""
    gauss_legendre_nodes(n, T=Float64)
    gauss_legendre01_nodes(n, T=Float64)

Return `n` Gauss-Legendre nodes on `[-1,1]` (or `[0,1]`), computed by Newton
iteration on `P_n`. Optimal interior quadrature nodes for polynomials of
degree up to `2n-1`.
"""
function gauss_legendre_nodes(n::Integer, ::Type{T}=Float64) where {T}
    x = [ cos(T(π) / (4n) * (4k - 2)) for k = n:-1:1 ]  # initial guess
    for _ = 1:100
        dx = legendre.(n, x) ./ dlegendre.(n, x)
        x .-= dx
        maximum(abs, dx) < 2 * eps(T) && return x
    end
    error("No convergence in Gauss-Legendre Newton iterations")
end

gauss_legendre01_nodes(n::Integer, ::Type{T}=Float64) where {T} = (gauss_legendre_nodes(n, T) .+ 1) ./ 2

"""
    gauss_legendre_quadrature(n, T=Float64)
    gauss_legendre01_quadrature(n, T=Float64)

Gauss-Legendre nodes and weights on `[-1,1]` (or `[0,1]`) with `n` points.
Degree of precision: `2n-1`.

```julia
x, w = gauss_legendre_quadrature(3)
w' * x.^4   # ≈ 2/5  (exact: DoP=5 ≥ 4)
```
"""
function gauss_legendre_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x = gauss_legendre_nodes(n, T)
    w = @. 2 / ((1 - x^2) * dlegendre(n, x)^2)
    x, w
end

function gauss_legendre01_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x, w = gauss_legendre_quadrature(n, T)
    (x .+ 1) ./ 2, w ./ 2
end

"""
    gauss_lobatto_nodes(n, T=Float64)
    gauss_lobatto01_nodes(n, T=Float64)

Return `n ≥ 2` Gauss-Lobatto nodes on `[-1,1]` (or `[0,1]`). Includes the
endpoints; interior nodes are computed by Newton iteration.
"""
function gauss_lobatto_nodes(n::Integer, ::Type{T}=Float64) where {T}
    n >= 2 || throw(ArgumentError("Gauss-Lobatto rules need at least 2 points"))
    n1 = n - 1
    x  = [ cos(T(π) * k / n1) for k = n1:-1:0 ]  # initial guess (Chebyshev)
    for _ = 1:100
        dx = @. (x * legendre(n1, x) - legendre(n1 - 1, x)) / ((n1 + 1) * legendre(n1, x))
        x .-= dx
        maximum(abs, dx) < 2 * eps(T) && return x
    end
    error("No convergence in Gauss-Lobatto Newton iterations")
end

gauss_lobatto01_nodes(n::Integer, ::Type{T}=Float64) where {T} = (gauss_lobatto_nodes(n, T) .+ 1) ./ 2

"""
    gauss_lobatto_quadrature(n, T=Float64)
    gauss_lobatto01_quadrature(n, T=Float64)

Gauss-Lobatto nodes and weights on `[-1,1]` (or `[0,1]`) with `n` points.
Degree of precision: `2n-3`.

```julia
x, w = gauss_lobatto_quadrature(4)
w' * x.^4   # ≈ 2/5  (exact: DoP=5 ≥ 4)
```
"""
function gauss_lobatto_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x = gauss_lobatto_nodes(n, T)
    w = @. 2 / ((n - 1) * n * legendre(n - 1, x)^2)
    x, w
end

function gauss_lobatto01_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x, w = gauss_lobatto_quadrature(n, T)
    (x .+ 1) ./ 2, w ./ 2
end

"""
    gauss_radau_nodes(n, T=Float64)
    gauss_radau01_nodes(n, T=Float64)

Return `n ≥ 2` left Gauss-Radau nodes on `[-1,1]` (or `[0,1]`): the left
endpoint plus `n-1` interior points, the roots of `P_n + P_{n-1}`. These are
the half-closed node sets used by some DG-SEM schemes; `1 .- reverse(x)`
gives the right-closed variant. They do not contain both vertices, so they
can serve as solver nodes but not as mesh nodes.
"""
function gauss_radau_nodes(n::Integer, ::Type{T}=Float64) where {T}
    n >= 2 || throw(ArgumentError("Gauss-Radau rules need at least 2 points"))
    x = [ -cos(2 * T(π) * k / (2n - 1)) for k = 0:n-1 ]  # initial guess; k = 0 is exactly -1
    for _ = 1:100
        f  = legendre.(n, x) .+ legendre.(n - 1, x)
        df = dlegendre.(n, x) .+ dlegendre.(n - 1, x)
        dx = f ./ df
        x .-= dx
        maximum(abs, dx) < 2 * eps(T) && return x
    end
    error("No convergence in Gauss-Radau Newton iterations")
end

gauss_radau01_nodes(n::Integer, ::Type{T}=Float64) where {T} = (gauss_radau_nodes(n, T) .+ 1) ./ 2

"""
    gauss_radau_quadrature(n, T=Float64)
    gauss_radau01_quadrature(n, T=Float64)

Left Gauss-Radau nodes and weights on `[-1,1]` (or `[0,1]`) with `n` points.
Degree of precision: `2n-2`.

```julia
x, w = gauss_radau_quadrature(3)
w' * x.^4   # ≈ 2/5  (exact: DoP=4)
```
"""
function gauss_radau_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x = gauss_radau_nodes(n, T)
    w = @. (1 - x) / (n^2 * legendre(n - 1, x)^2)
    x, w
end

function gauss_radau01_quadrature(n::Integer, ::Type{T}=Float64) where {T}
    x, w = gauss_radau_quadrature(n, T)
    (x .+ 1) ./ 2, w ./ 2
end

###########################################################################
## Multi-dimensional quadrature

"""
    quadrature(eg::ElementGeometry, p, T=Float64)

Return `(ξ, w)`: quadrature nodes `ξ` (`npts × D`) and weights `w` (`npts`)
exact for polynomials of degree `p` on the reference element `eg`, with number
type `T`.

- `Block`: tensor-product Gauss-Legendre rule on `[0,1]^D` with `cld(p+1, 2)`
  points per direction.
- `Simplex`: symmetric rules tabulated in `simplex_quadrature_rules.txt`,
  looked up by `(D, p)`. See [`max_quadrature_degree`](@ref).

The weights sum to the volume of the reference element: 1 for blocks and
`1/factorial(D)` for simplices.
"""
function quadrature(eg::Block{D}, p::Integer, ::Type{T}=Float64) where {D,T}
    n = cld(p + 1, 2)
    ξ0, w0 = gauss_legendre01_quadrature(n, T)
    ξ = tensor_nodes(eg, ξ0)
    w = [ prod(w0[i[d]] for d in 1:D; init=one(T))
          for i in vec(collect(Iterators.product(ntuple(_ -> 1:n, D)...))) ]
    ξ, w
end

const _simplex_rules = Dict{Tuple{Int,Int},Tuple{Matrix{Float64},Vector{Float64}}}()

# Parse src/basis/simplex_quadrature_rules.txt into _simplex_rules, once.
function _load_simplex_rules!()
    file = joinpath(@__DIR__, "simplex_quadrature_rules.txt")
    open(file) do io
        for line in eachline(io)
            startswith(line, "#") && continue
            isempty(strip(line)) && continue
            tag, D, p, n = split(line)
            tag == "rule" || error("Malformed simplex quadrature rule file: $file")
            D, p, n = parse(Int, D), parse(Int, p), parse(Int, n)
            ξ = zeros(n, D)
            w = zeros(n)
            for i in 1:n
                vals = parse.(Float64, split(readline(io)))
                ξ[i, :] = vals[1:D]
                w[i] = vals[D+1]
            end
            _simplex_rules[(D, p)] = (ξ, w)
        end
    end
end

"""
    max_quadrature_degree(::Simplex{D}) where {D}

Highest polynomial degree for which a tabulated simplex quadrature rule
exists in dimension `D`.
"""
function max_quadrature_degree(::Simplex{D}) where {D}
    isempty(_simplex_rules) && _load_simplex_rules!()
    maximum(p for (d, p) in keys(_simplex_rules) if d == D)
end

function quadrature(::Simplex{D}, p::Integer, ::Type{T}=Float64) where {D,T}
    isempty(_simplex_rules) && _load_simplex_rules!()
    haskey(_simplex_rules, (D, p)) ||
        error("No simplex quadrature rule for D=$D, p=$p; maximum degree is $(max_quadrature_degree(Simplex{D}()))")
    ξ, w = _simplex_rules[(D, p)]
    T.(ξ), T.(w)
end
