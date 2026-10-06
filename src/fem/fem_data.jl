###########################################################################
## Precomputed data for finite element assembly
#
# Two kinds of data, separated by lifetime and size:
#
#   RefOps      the reference element evaluated at a point set (the volume
#               quadrature points, or the quadrature points of one face):
#               tiny and mesh independent.
#   Metric      metric terms at the quadrature points of every element:
#               w·det J, J⁻¹ and x. A few times the size of the mesh.
#   FaceMetric  the same at the face quadrature points of every face.
#
# No physical basis gradients are stored: they are applied on the fly in the
# kernels from ϕξ and J⁻¹, which is as fast and 30-300 times smaller.
#
# Layout: reference gradients are stored as an (ng·D) × ns matrix ϕξ whose
# row g + ng(a-1) holds ∂ϕ_i/∂ξ_a at point g, so that a sum over points and
# directions is one matrix product. The D-vector of one point is read and
# written with gradvec / setgradvec!. Coordinates follow the mesh
# convention J[i,a] = ∂x_i/∂ξ_a, so a physical gradient is ∇ₓϕ = ∇ξϕ·J⁻¹.

"""
    RefOps{D,T}

Reference element evaluated at a set of `ng` points with weights:

- `w`:  quadrature weights (`ng`)
- `ϕ`:  shape functions at the points (`ng × ns`)
- `ϕξ`: reference gradients as an `(ng·D) × ns` matrix; row `g + ng(a-1)` holds
  `∂ϕ_i/∂ξ_a` at point `g`, so that a sum over points and directions is one
  matrix product.

`RefOps(fe, qdeg)` uses the volume quadrature rule of degree `qdeg`;
`RefOps(fe, ξ, w)` an explicit point set. See [`face_refops`](@ref) for the
trace operators on the faces.
"""
struct RefOps{D,T}
    w::Vector{T}
    ϕ::Matrix{T}
    ϕξ::Matrix{T}
end

"""Number of points of a `RefOps`."""
npoints(r::RefOps) = length(r.w)

function RefOps(fe::FiniteElement{D,G,T}, ξ::AbstractMatrix, w::AbstractVector) where {D,G,T}
    dϕ = dshapefcns(fe, ξ)                                     # ng × ns × D
    RefOps{D,T}(T.(w), shapefcns(fe, ξ), reshape(permutedims(dϕ, (1,3,2)), length(w)*D, :))
end

RefOps(fe::FiniteElement{D,G}, qdeg::Integer) where {D,G} = RefOps(fe, quadrature(G(), qdeg)...)

"""
    face_refops(fe, qdeg) -> Vector{RefOps}

Trace operators on each local face of `fe`: the face quadrature rule of
degree `qdeg` mapped into the volume reference coordinates. The weights are
those of the face reference element; the physical scaling is in
[`FaceMetric`](@ref).
"""
function face_refops(fe::FiniteElement{D,G,T}, qdeg::Integer) where {D,G,T}
    eg = G()
    s, w = quadrature(subgeom(eg, D-1), qdeg)
    N  = shapefcns(subgeom(eg, D-1), s)                        # linear face map (ngf × nfv)
    V  = vertices(eg)
    [ RefOps(fe, N * V[fm, :], w) for fm in eachcol(facemap(eg)) ]
end

# The D-vector stored at rows g, g+ng, ..., g+(D-1)ng of column i
@inline gradvec(A, g, i, ng, ::Val{D}) where {D} =
    SVector{D}(ntuple(a -> @inbounds(A[g + ng*(a-1), i]), Val(D)))
@inline function setgradvec!(A, v::SVector{D}, g, i, ng) where {D}
    @inbounds for a in 1:D
        A[g + ng*(a-1), i] = v[a]
    end
end

###########################################################################
## Metric terms

"""
    jacobians(m, r::RefOps) -> (J, x)

Jacobian `J[g,e]` (an `SMatrix`, `J[i,a] = ∂x_i/∂ξ_a`) and physical coordinates
`x[g,e]` (an `SVector`) at the points of `r` for every element of `m`, as two
`ng × nel` arrays.
"""
function jacobians(m::HighOrderMesh{D,G,T}, r::RefOps{D,T}) where {D,G,T}
    ng, ns, ne = npoints(r), nnodes(m.fe), nel(m)
    X  = reshape(dg_nodes(m), ns, ne*D)                        # column e + ne(i-1): x_i of element e
    dX = r.ϕξ * X                                              # (ng·D) × (ne·D)
    xq = r.ϕ  * X                                              # ng × (ne·D)
    J  = [ SMatrix{D,D,T}(ntuple(k -> dX[g + ng*((k-1)÷D), e + ne*((k-1)%D)], D*D)) for g in 1:ng, e in 1:ne ]
    x  = [ SVector{D,T}(ntuple(i -> xq[g, e + ne*(i-1)], D)) for g in 1:ng, e in 1:ne ]
    J, x
end

"""
    Metric{D,T,L}

Metric terms at the quadrature points of every element, each an `ng × nel`
array: `wJ` is `w·det J`, `Jinv` the inverse Jacobian (`SMatrix`) and `x` the
physical coordinates (`SVector`). `Metric(m, r::RefOps)` builds it.
"""
struct Metric{D,T,L}
    wJ::Matrix{T}
    Jinv::Matrix{SMatrix{D,D,T,L}}
    x::Matrix{SVector{D,T}}
end

function Metric(m::HighOrderMesh{D,G,T}, r::RefOps{D,T}) where {D,G,T}
    J, x = jacobians(m, r)
    detJ = det.(J)
    all(>(0), detJ) || error("the mesh has elements with a non-positive Jacobian")
    Metric(r.w .* detJ, inv.(J), x)
end

# Outward normal of local face j of the reference element, scaled so that
# w·dS·n = w_face · det J · J⁻ᵀ ν  (Nanson's formula) with w_face the weights
# of the face reference element.
function refnormal(eg::ElementGeometry{D}, j) where {D}
    V  = vertices(eg)
    v  = [ SVector{D,Float64}(V[k, :]) for k in facemap(eg)[:, j] ]
    ν  = D == 1 ? SVector(1.0) :
         D == 2 ? (t = v[2] - v[1]; SVector(t[2], -t[1])) :
                  cross(v[2] - v[1], v[3] - v[1])
    outward = dot(ν, sum(v)/length(v) - SVector{D}(sum(V, dims=1)/size(V,1)))
    outward > 0 ? ν : -ν
end

"""
    FaceMetric{D,T,L}

Metric terms at the face quadrature points of every face of every element,
each an `ngf × nfaces × nel` array: `wn` is `w·dS·n` with `n` the unit outward
normal (`SVector`), `Jinv` the inverse Jacobian and `x` the physical
coordinates. `FaceMetric(m, fr)` builds it from the trace operators `fr` of
[`face_refops`](@ref).
"""
struct FaceMetric{D,T,L}
    wn::Array{SVector{D,T},3}
    Jinv::Array{SMatrix{D,D,T,L},3}
    x::Array{SVector{D,T},3}
end

function FaceMetric(m::HighOrderMesh{D,G,T}, fr::Vector{RefOps{D,T}}) where {D,G,T}
    nf, ne = size(m.nb)
    ngf    = npoints(fr[1])
    wn   = Array{SVector{D,T},3}(undef, ngf, nf, ne)
    Jinv = Array{SMatrix{D,D,T,D*D},3}(undef, ngf, nf, ne)
    x    = Array{SVector{D,T},3}(undef, ngf, nf, ne)
    for j in 1:nf
        J, xj = jacobians(m, fr[j])
        ν = refnormal(G(), j)
        wn[:, j, :]   = [ fr[j].w[g] * det(J[g,e]) * (inv(J[g,e])' * ν) for g in 1:ngf, e in 1:ne ]
        Jinv[:, j, :] = inv.(J)
        x[:, j, :]    = xj
    end
    FaceMetric(wn, Jinv, x)
end

# nbpt[g,j,e]: index on the neighbor's face of the point coinciding with point g
# of face j of element e (0 on boundary faces). Points are matched by position
# relative to the face centroid, so that periodic faces match too.
function neighbor_points(m::HighOrderMesh, fm::FaceMetric)
    ngf, nf, ne = size(fm.x)
    nbpt = zeros(Int, ngf, nf, ne)
    for e in 1:ne, j in 1:nf
        nb = m.nb[j, e]
        isboundary(nb) && continue
        x1 = view(fm.x, :, j, e)
        x2 = view(fm.x, :, nb.face, nb.el)
        c  = (sum(x1) - sum(x2)) / ngf
        for g in 1:ngf
            dist = [ norm(x1[g] - c - x2[h]) for h in 1:ngf ]
            k = argmin(dist)
            dist[k] <= 1e-8 * (1 + norm(x1[g])) || error("face points of elements $e and $(nb.el) do not coincide")
            nbpt[g, j, e] = k
        end
    end
    nbpt
end

###########################################################################
## Assembly data

"""
    FEMData{D}

Abstract supertype of [`CGData`](@ref) and [`DGData`](@ref). Both hold the
mesh `m`, the volume reference operators `ref` and the metric terms `met`, on
which the elemental kernels of `assembly.jl` dispatch.
"""
abstract type FEMData{D} end

"""
    CGData(m; qdeg=2p+1)

Volume data for continuous Galerkin assembly on the mesh `m`: the reference
element at the quadrature points of degree `qdeg` (`ref::RefOps`) and the
metric terms of every element (`met::Metric`).
"""
struct CGData{D,G,T,L} <: FEMData{D}
    m::HighOrderMesh{D,G,T}
    ref::RefOps{D,T}
    met::Metric{D,T,L}
end

function CGData(m::HighOrderMesh; qdeg=2porder(m)+1)
    ref = RefOps(m.fe, qdeg)
    CGData(m, ref, Metric(m, ref))
end

"""
    DGData(m; qdeg=2p+1)

Data for discontinuous Galerkin assembly on the mesh `m`: the volume data of
[`CGData`](@ref) plus, for every face of every element, the trace operators
(`fref[j]` for local face `j`), the face metric terms (`fmet::FaceMetric`)
and the matching of the face points with the neighbor's (`nbpt[g,j,e]`, the
index on the neighbor's face of the point coinciding with point `g` of face
`j` of element `e`).
"""
struct DGData{D,G,T,L} <: FEMData{D}
    m::HighOrderMesh{D,G,T}
    ref::RefOps{D,T}
    met::Metric{D,T,L}
    fref::Vector{RefOps{D,T}}
    fmet::FaceMetric{D,T,L}
    nbpt::Array{Int,3}
end

function DGData(m::HighOrderMesh; qdeg=2porder(m)+1)
    ref  = RefOps(m.fe, qdeg)
    fref = face_refops(m.fe, qdeg)
    fmet = FaceMetric(m, fref)
    DGData(m, ref, Metric(m, ref), fref, fmet, neighbor_points(m, fmet))
end

Base.show(io::IO, d::CGData) =
    print(io, "CGData: $(nel(d.m)) elements, $(nnodes(d.m.fe)) nodes and $(npoints(d.ref)) quadrature points per element")
Base.show(io::IO, d::DGData) =
    print(io, "DGData: $(nel(d.m)) elements, $(nnodes(d.m.fe)) nodes, $(npoints(d.ref)) volume and $(npoints(d.fref[1])) face quadrature points per element")
