###########################################################################
## Polynomials

function legendre_poly(x::AbstractVector{T}, p::Int; diff=false) where T
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

function multivar_legendre_poly(x::AbstractArray{T}, p::Int) where T
    D = size(x,2)
    ix = [ i for i in Iterators.product(fill(1:p+1,D)...)][:]
    Ls = [ legendre_poly(xx, p) for xx in eachcol(x) ]

    P = zeros(T, size(x,1), length(ix))
    for i = 1:size(x,1)
        for j = 1:length(ix)
            val = one(T)
            for k = 1:D
                val *= Ls[k][i,ix[j][k]]
            end
            P[i,j] = val
        end
    end
    P
end

function multivar_monomial_poly(x::AbstractArray{T}, p::Int) where T
    isempty(x) && return [one(T);;]
    x = x'
    D = size(x,1)
    powers = sort([ i for i in Iterators.product(fill(0:p,D)...) if sum(i) <= p ],
                  lt=(i,j)->sum(i)<sum(j))
    [ prod(xx .^ k) for xx in eachcol(x), k in powers ]
end

eval_poly(eg::Simplex, s, p) = multivar_monomial_poly(s, p)
eval_poly(eg::Block, s, p) = multivar_legendre_poly(s, p)

###########################################################################
## Nodes

equispaced(p::Int) = (0:p) // p

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

eval_poly(::FiniteElement{D,G,P,T}, s) where {D,G,P,T} = eval_poly(G(), s, P)

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

eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}) where {D,G,P,T} =
    eval_poly(G(), ss, P) * fe.shapefcn_coeff.vol

eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}) where {T} =
    eval_shapefcns(FiniteElement(eg, 1, T), ss)

function eval_fcn(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}) where {D,G,P,T}
    nss,ndim = size(ss,1),size(ss,2)
    ns,nel = size(u,1),size(u,2)
    geomdim = ndim2geomdim(D, ndim)
    C = fe.shapefcn_coeff[geomdim]

    matrix_output = eval_poly(fe, ss) * (C * reshape(u,ns,:))
    if ndims(u) <= 2
        return matrix_output
    else
        return reshape(matrix_output, nss, nel, :)
    end
end

#################################################################################
