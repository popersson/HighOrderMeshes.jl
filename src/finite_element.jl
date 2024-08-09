###########################################################################
## FiniteElement

struct FiniteElement{D,G<:ElementGeometry{D},P,T}
    ref_nodes::NTuple{D, Matrix{T}}
    shapefcn_coeff::NTuple{D, Matrix{T}}
end

ref_nodes(fe::FiniteElement{D,G,P,T}, dim) where {D,G,P,T} =
    dim == 0 ? zeros(T,1,0) : fe.ref_nodes[dim]

"""
    eval_poly(::FiniteElement{D,G,P,T}, s; gradient=false) where {D,G,P,T}

TBW
"""
eval_poly(::FiniteElement{D,G,P,T}, s; gradient=false) where {D,G,P,T} = eval_poly(G(), s, P; gradient=gradient)

elgeom(::FiniteElement{D,G,P,T}) where {D,G,P,T} = G()
dim(::FiniteElement{D,G,P,T}) where {D,G,P,T} = D
porder(::FiniteElement{D,G,P,T}) where {D,G,P,T} = P
nbr_ho_nodes(fe::FiniteElement{D}, dims=D) where {D} = size(ref_nodes(fe,dims),1)
function corner_nodes(fe::FiniteElement{D,G,P,T}) where {D,G,P,T}
    fe1 = FiniteElement(G(), 1)
    indexin(eachrow(ref_nodes(fe1,D)), eachrow(ref_nodes(fe,D)))
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
    nodes = Tuple( ref_nodes(subgeom(eg,d), sline) for d in 1:D )
    coeffs = Tuple( inv(eval_poly(eg, nodes[d], p)) for d in 1:D )
    FiniteElement{D,typeof(eg),p,T}(nodes, coeffs)
end

FiniteElement(eg::ElementGeometry, p::Int, T=Float64) = FiniteElement(eg, T.(equispaced(p)))

###########################################################################
## Shape functions

"""
    eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}; gradient=false) where {T}

    Evaluate shape functions at reference coordinates ss

    Dimensions:
      ss: nss x ndim
      output, gradient=false: nss x ns
      output, gradient=true:  nss x ns x ndim
"""
function eval_shapefcns(fe::FiniteElement{D,G,P,T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
    nss,ndim = size(ss,1),size(ss,2)
    ndim == 0 && return ones(T,size(ss,1))
    pol = eval_poly(G(), ss, P; gradient=gradient)
    C = fe.shapefcn_coeff[ndim]
    if gradient
        return cat( (pol[:,:,k] * C for k = 1:ndim)..., dims=3)
    else
        return pol * C
    end
end    

eval_shapefcns(eg::ElementGeometry, ss::AbstractArray{T}; gradient=false) where {T} =
    eval_shapefcns(FiniteElement(eg, 1, T), ss; gradient=gradient)

"""
    eval_fcn(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}

    Evaluate finite element function fe,u at reference coordinates ss

    Dimensions:
      u: ns x nel x nc  or  ns x nel
      ss: nss x ndim
      output, gradient=false: nss x nel x nc  (3rd dim dropped if u is ns x nel)
      output, gradient=true:  nss x nel x nc x ndim
"""
function eval_fcn(fe::FiniteElement{D,G,P,T}, u::Array{T}, ss::AbstractArray{T}; gradient=false) where {D,G,P,T}
   
    nss,ndim = size(ss,1),size(ss,2)
    ns,nel = size(u,1),size(u,2)
    @assert ns == nbr_ho_nodes(fe, ndim)

    ndim == 0 && return gradient ? zeros(T,nss,nel,nc,0) : repeat(u[[1],:,:],nss)
    
    C = fe.shapefcn_coeff[ndim]

    V = eval_poly(fe, ss; gradient=gradient)
    if gradient
        V = reshape(permutedims(V,(1,3,2)),:,ns)
    end
    matrix_output = V * (C * reshape(u,ns,:))
    if gradient
        return permutedims(reshape(matrix_output, nss, ndim, nel, size(u,3)), (1,3,4,2))
    else
        if ndims(u) <= 2
            return matrix_output
        else
            return reshape(matrix_output, nss, nel, :)
        end
    end
end

#################################################################################
