using LinearAlgebra: eigen, SymTridiagonal

"""
    gauss_legendre_quadrature(p, T::Type=Float64)

Gauss-Legendre quadrature on [-1,1], degree of precision p.

Example usage:
```julia
# Integrate x^4 on [-1,1]
x,w = gauss_legendre_quadrature(4)
w' * x.^4
```
"""
function gauss_legendre_quadrature(p::Int; T::Type=Float64)
    n = p÷2 + 1
    b = 1:n-1
    b = @. b / sqrt(T(4*b^2 - 1))
    eval, evec = eigen(SymTridiagonal(zeros(T,n), b))
    return (eval .+ 1) ./ 2, evec[1,:].^2
end

function gauss_lobatto_nodes(n1::Int; T::Type=Float64)
    n = n1 - 1
    x = @. cos(pi*((0:n)/T(n)))
    for it = 1:100
        P = legendre_poly_neg1pos1(x,n)
        xold = copy(x)
        x = @. xold - (x * P[:,end] - P[:,end-1]) / ( (n+1)*P[:,end] )
        if maximum(abs.(x - xold)) < 2*eps(T)
            return (reverse(x) .+ 1) ./ 2
        end
    end
    throw("No convergence in Gauss-Lobatto Newton iterations")
end

function gauss_lobatto_quadrature(p::Int; T::Type=Float64)
    n = (p+2)÷2 + 1
    x = gauss_lobatto_nodes(n, T=T)
    P = legendre_poly(x,n)
    w = @. 1 / (n * (n+1) * P[:,end]^2)
    return x, w
end

quadrature(eg::ElementGeometry, p) = throw("Quadrature not implemented")

quadrature(eg::Simplex{D}, p) where D = simplex_quadrature(SimplexQuadRule{D,p}())

function quadrature(eg::Block{D}, p; T::Type=Float64) where D
    ξ0,w0 = gauss_legendre_quadrature(p, T=T)
    ξ = ref_nodes(eg, ξ0)
    w = w0
    for d = 2:D
        w = w .* w0'
        w = w[:]
    end
    ξ,w
end
