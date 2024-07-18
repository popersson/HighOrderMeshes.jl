###########################################################################
## Misc utils

"""
    matrix_tensor_product(A,B)

TBW
"""
matrix_tensor_product(A,B) = reshape(A * reshape(B, size(A,2), :), size(A,1), size(B)[2:end]...)

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
    if gradient
        dLs = [ legendre_poly(xx, p, diff=true)[2] for xx in eachcol(x) ]
        P = zeros(T, size(x,1), length(ix), D)
        for k = 1:D
            cLs = [ i == k ? dLs[i] : Ls[i] for i = 1:D ]
            make_outer_prod!(cLs, view(P, :, :, k))
        end
    else
        P = zeros(T, size(x,1), length(ix))
        make_outer_prod!(Ls, P)
    end

    P
end

"""
    multivar_monomial_poly(x::AbstractArray{T}, p::Int; gradient=false) where T

TBW
"""
function multivar_monomial_poly(x::AbstractArray{T}, p::Int; gradient=false) where T
    ### TODO: Gradient
    isempty(x) && return [one(T);;]
    x = x'
    D = size(x,1)
    powers = sort([ i for i in Iterators.product(fill(0:p,D)...) if sum(i) <= p ],
                  lt=(i,j)->sum(i)<sum(j))
    [ prod(xx .^ k) for xx in eachcol(x), k in powers ]
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
    ref_nodes(eg::Simplex{D}, s1) where {D}

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

function ref_nodes(eg::Simplex{D}, s1) where {D}
    sD = ref_nodes(Block{D}(), s1)
    sD = sD[sum(sD,dims=2)[:,1].<=1,:]
end

###########################################################################
## FiniteElement

struct FiniteElement{D,G<:ElementGeometry{D},P,T}
    ref_nodes::@NamedTuple{line::AbstractVector{T},
                           face::AbstractArray{T},
                           vol::AbstractArray{T}}
    shapefcn_coeff::@NamedTuple{line::AbstractArray{T},
                                face::AbstractArray{T},
                                vol::AbstractArray{T}}
end

ndim2geomdim(D, ndim) = D == ndim ? :vol : D-1 == ndim ? :face : 1 == ndim : :line : throw("Geometry dimension not implemented")

"""
    eval_poly(::FiniteElement{D,G,P,T}, s; gradient=false) where {D,G,P,T}

TBW
"""
eval_poly(::FiniteElement{D,G,P,T}, s; gradient=false) where {D,G,P,T} = eval_poly(G(), s, P; gradient=gradient)

elgeom(::FiniteElement{D,G,P,T}) where {D,G,P,T} = G()
dim(::FiniteElement{D,G,P,T}) where {D,G,P,T} = D
porder(::FiniteElement{D,G,P,T}) where {D,G,P,T} = P
nbr_ho_nodes(fe::FiniteElement) = size(fe.ref_nodes.vol,1)
function corner_nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T}
    fe1 = FiniteElement(G(), 1)
    indexin(eachrow(fe1.ref_nodes.vol), eachrow(fe.ref_nodes.vol))
end

name(::FiniteElement{D,G,P,T}) where {D,G,P,T} = "p=" * string(P) * " " * name(G())

function Base.show(io::IO, fe::FiniteElement)
    print(io, "FiniteElement: $(dim(fe))D, $(name(fe)) element")
end

"""
    FiniteElement(eg::ElementGeometry{D}, sline::AbstractVector{T}) where {D,T}
    FiniteElement(eg::ElementGeometry, p::Int, T=Float64)

TBW
"""
function FiniteElement(eg::ElementGeometry{D}, sline::AbstractVector{T}) where {D,T}
    p = length(sline) - 1
    sface = ref_nodes(facegeom(eg), sline)
    svol = ref_nodes(eg, sline)

    Cline = inv(eval_poly(eg, sline, p))
    Cface = inv(eval_poly(eg, sface, p))
    Cvol = inv(eval_poly(eg, svol, p))

    FiniteElement{D,typeof(eg),p,T}((line=sline,face=sface,vol=svol),
                                    (line=Cline,face=Cface,vol=Cvol))
end

FiniteElement(eg::ElementGeometry, p::Int, T=Float64) = FiniteElement(eg, T.(equispaced(p)))

###########################################################################
## Shape functions

"""
    eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}; gradient=false) where {T}

TBW
"""
function eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    pol = eval_poly(G(), ss, P; gradient=gradient)
    C = fe.shapefcn_coeff[:vol]
    if gradient
        return cat( (pol[:,:,k] * C for k = 1:D)..., dims=3)
    else
        return pol * C
    end
end    

eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}; gradient=false) where {T} =
    eval_shapefcns(FiniteElement(eg, 1, T), ss; gradient=gradient)

"""
    eval_fcn(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}

TBW
"""
function eval_fcn(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    nss,ndim = size(ss,1),size(ss,2)
    ns,nel = size(u,1),size(u,2)
    geomdim = ndim2geomdim(D, ndim)
    C = fe.shapefcn_coeff[geomdim]

    V = eval_poly(fe, ss; gradient=gradient)
    if gradient
        V = reshape(permutedims(V,(1,3,2)),:,ns)
    end
    matrix_output = V * (C * reshape(u,ns,:))
    if gradient
        return permutedims(reshape(matrix_output, nss, D, nel, size(u,3)), (1,3,4,2))
    else
        if ndims(u) <= 2
            return matrix_output
        else
            return reshape(matrix_output, nss, nel, :)
        end
    end
end

#################################################################################
