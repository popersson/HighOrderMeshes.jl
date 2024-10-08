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
