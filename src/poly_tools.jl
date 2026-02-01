###########################################################################
## Polynomials

"""
    legendre_poly_neg1pos1(x::AbstractVector{T}, p::Int; diff=false) where T

TBW
"""
function legendre_poly_neg1pos1(x::AbstractVector{T}, p::Int; diff=false) where T
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
    legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T

TBW
"""
function legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T
    # Legendre on [0,1]
    ys = legendre_poly_neg1pos1(2x .- 1, p; diff=diff)
    if diff
        ys[2] .*= 2
    end
    ys
end

"""
    multivar_legendre_poly(x::AbstractArray{T}, p::Int; gradient=false) where T

TBW
"""
function multivar_legendre_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
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

    Ls = [ legendre_poly(xx, p) for xx in eachcol(x) ]
    if !gradient
        P = zeros(T, size(x,1), length(ix))
        make_outer_prod!(Ls, P)
    else
        dLs = [ legendre_poly(xx, p, diff=true)[2] for xx in eachcol(x) ]
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

TBW
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
eval_poly(eg::Block, s, p; gradient=false) = multivar_legendre_poly(s, p; gradient=gradient)

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
