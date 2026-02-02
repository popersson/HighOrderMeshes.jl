###########################################################################
## Polynomials

"""
    legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T

Legendre polynomials L_p(x) on [-1,1], vectorized over evaluation points in x.
Optional output derivative L'_p(x).
"""
function legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T
    # Legendre on [-1,1]
    y = ones(T, length(x), p+1)
    y[:,2] .= x
    if diff
        dy = zeros(T, length(x), p+1)
        dy[:,2] .= one(T)
    end
    for j = 1:p-1
        @. @views y[:,j+2] = ((2j+1)*x*y[:,j+1] - j*y[:,j]) / (j+1)
        if diff
            @. @views dy[:,j+2] = ((2j+1)*(x*dy[:,j+1] + y[:,j+1]) - j*dy[:,j]) / (j+1)
        end
    end
    diff ? (y,dy) : y
end

"""
    legendre01_poly(x::AbstractVector{T}, p::Int; diff=false) where T

Shifted Legendre polynomials on [0,1].
"""
function legendre01_poly(x::AbstractVector{T}, p::Int; diff=false) where T
    # Legendre on [0,1]
    ys = legendre_poly(2x .- 1, p; diff=diff)
    if diff
        ys[2] .*= 2
    end
    ys
end

"""
    multivar_legendre01_poly(x::AbstractArray{T}, p::Int; gradient=false) where T

Multivariate Shifted Legendre polynomials of degree p on [0,1]^D where D = size(x,2).
"""
function multivar_legendre01_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
    D = size(x,2)
    ix = [ i for i in Iterators.product(fill(1:p+1,D)...)][:]
    function make_outer_prod!(Ls, P)
        for i = 1:size(x,1)
            for j = 1:length(ix)
                val = one(T)
                for k = 1:D
                    val *= Ls[k][i,ix[j][k]]
                end
                P[i,j] = val
            end
        end
    end

    Ls = [ legendre01_poly(xx, p) for xx in eachcol(x) ]
    if !gradient
        P = zeros(T, size(x,1), length(ix))
        make_outer_prod!(Ls, P)
    else
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
    multivar_monomial_poly(x::AbstractArray{T}, p::Int; gradient=false) where T

Multivariate monomials of degree p on R^D where D = size(x,2).
"""
function multivar_monomial_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
    isempty(x) && return [one(T);;]
    x = x'
    D = size(x,1)
    powers = sort([ i for i in Iterators.product(fill(0:p,D)...) if sum(i) <= p ],
                  lt=(i,j)->sum(i)<sum(j))
    if !gradient
        return [ prod(xx .^ k) for xx in eachcol(x), k in powers ]
    else
        p = zeros(T, size(x,2), length(powers), D)
        for d = 1:D
            mul = (ix->ix[d]).(powers)
            mask = Tuple( Int(i==d) for i = 1:D)
            diffpowers = [ max.(0, power .- mask) for power in powers ]
            p[:,:,d] = [ m .* prod(xx .^ k) for xx in eachcol(x), (k,m) in zip(diffpowers,mul) ]
        end
        p
    end
end

"""
    eval_poly(eg::Simplex, s, p; gradient=false)
    eval_poly(eg::Block, s, p; gradient=false)

TBW
"""
eval_poly(eg::Simplex, s, p; gradient=false) = multivar_monomial_poly(s, p; gradient=gradient)
eval_poly(eg::Block, s, p; gradient=false) = multivar_legendre01_poly(s, p; gradient=gradient)

###########################################################################
## 1D Quadrature

function gauss_legendre_nodes(n::Int; T::Type=Float64)
    x = [ cos(T(π) / 4n * (4k - 2)) for k = n:-1:1 ]
    for it = 1:100
        P,dP = legendre_poly(x, n, diff=true)
        xold = copy(x)
        @. x = xold - P[:,end] / dP[:,end]
        if maximum(abs.(x - xold)) < 2*eps(T)
            return x
        end
    end
    throw("No convergence in Gauss-Legendre Newton iterations")
end

gauss_legendre01_nodes(n::Int; T::Type=Float64) = ( gauss_legendre_nodes(n, T=T) .+ 1) ./ 2

"""
    gauss_legendre_quadrature(n; T::Type=Float64)

Gauss-Legendre quadrature on [-1,1] with n points.
Degree of Precision = 2n - 1

Example usage:
```julia
# Compute integral(x^4, x=[-1,1]) = 2/5
# n = 3 => DoP = 5 => Exact
x,w = gauss_legendre_quadrature(3)
w' * x.^4
```
"""
function gauss_legendre_quadrature(n::Int; T::Type=Float64)
    # Degree of precision = 2n - 1
    x = gauss_legendre_nodes(n, T=T)
    P,dP = legendre_poly(x, n, diff=true)
    w = @. 2 / (1 - x^2) / dP[:,end]^2
    return x,w
end

function gauss_legendre01_quadrature(n::Int; T::Type=Float64)
    x,w = gauss_legendre_quadrature(n, T=T)
    (x .+ 1) ./ 2, w ./ 2
end

function gauss_lobatto_nodes(n::Int; T::Type=Float64)
    n1 = n - 1
    x = @. cos(pi*((n1:-1:0)/T(n1)))
    for it = 1:100
        P = legendre_poly(x,n1)
        xold = copy(x)
        x = @. xold - (x * P[:,end] - P[:,end-1]) / ( (n1+1)*P[:,end] )
        if maximum(abs.(x - xold)) < 2*eps(T)
            return x
        end
    end
    throw("No convergence in Gauss-Lobatto Newton iterations")
end

gauss_lobatto01_nodes(n::Int; T::Type=Float64) = ( gauss_lobatto_nodes(n, T=T) .+ 1) ./ 2

"""
    gauss_lobatto_quadrature(n; T::Type=Float64)

Gauss-Lobatto quadrature on [-1,1] with n points.
Degree of Precision = 2n - 3

Example usage:
```julia
# Compute integral(x^4, x=[-1,1]) = 2/5
# n = 4 => DoP = 5 => Exact
x,w = gauss_lobatto_quadrature(4)
w' * x.^4
```
"""
function gauss_lobatto_quadrature(n::Int; T::Type=Float64)
    # Degree of precision = 2n - 3
    x = gauss_lobatto_nodes(n, T=T)
    P = legendre_poly(x, n-1)
    w = @. 2 / ((n-1) * n * P[:,end]^2)
    return x, w
end

function gauss_lobatto01_quadrature(n::Int; T::Type=Float64)
    x,w = gauss_lobatto_quadrature(n, T=T)
    (x .+ 1) ./ 2, w ./ 2
end

###########################################################################
## Nodes

equispaced(p::Int) = (0:p) // p

"""
    ref_nodes(::Block{D}, s1::AbstractVector{T}) where {D,T}
    ref_nodes(::Simplex{D}, s1) where {D}

TBW
"""
function ref_nodes(::Block{D}, s1::AbstractVector{T}) where {D,T}
    sD = zeros(T, 1, 0)
    for d = 1:D
        sD = hcat(repeat(sD,length(s1),1),
                  reshape(repeat(s1',size(sD,1),1),:,1))
    end
    sD
end

function ref_nodes(::Simplex{D}, s1) where {D}
    sD = ref_nodes(Block{D}(), s1)
    sD = sD[sum(sD,dims=2)[:,1].<=1,:]
end

## Identify geometry from number of high-order nodes
function find_elgeom(D, P, nnodes)
    if (P+1)^D == nnodes
        return Block{D}()
    elseif binomial(P+D,D) == nnodes
        return Simplex{D}()
    else
        error("Cannot determine element shape")
    end
end

###########################################################################
## General Quadrature

quadrature(eg::ElementGeometry, p) = throw("Quadrature not implemented")

quadrature(eg::Simplex{D}, p) where D = simplex_quadrature(SimplexQuadRule{D,p}())

function quadrature(eg::Block{D}, p; T::Type=Float64) where D
    ξ0,w0 = gauss_legendre01_quadrature(p, T=T)
    
    ξ = ref_nodes(eg, ξ0)
    w = w0
    for d = 2:D
        w = w .* w0'
        w = w[:]
    end
    ξ,w
end

