using Test

# -------------------------------------------------------------------
# 1. Core Library Tests
# -------------------------------------------------------------------
module TestCore
using Test
using HighOrderMeshes
using LinearAlgebra
import DistMesh
const HOM = HighOrderMeshes

# Evaluate the monomial ξ^α at every row of ξ.
monomial(ξ, α) = [ prod(ξ[k,d]^α[d] for d in eachindex(α); init=1.0) for k in axes(ξ,1) ]

# Volume of a mesh by quadrature of the Jacobian determinant, and the smallest
# ratio of the determinant to its maximum in the same element (> 0: not inverted).
function volume_and_quality(m)
    ξ, w = quadrature(elgeom(m), 3*porder(m) + 2)
    J = interpolate(dshapefcns(m.fe, ξ), dg_nodes(m))
    dets = [ det(J[g,e,:,:]) for g in axes(J,1), e in axes(J,2) ]
    sum(w .* dets), minimum(dets ./ maximum(dets, dims=1))
end

# Coordinates of the nodes of every boundary face with tag `tag`, one matrix per face
boundary_face_nodes(m, tag) = [ m.x[m.el[mkface2nodes(m)[:,j],iel],:]
    for j in axes(m.nb,1), iel in axes(m.nb,2) if isboundary(m.nb[j,iel]) && bndtag(m.nb[j,iel]) == tag ]

# Central finite-difference gradient of f(ξ) (f returns a matrix) along coordinate d.
function fdgrad(f, ξ, d; h=1e-6)
    e = zeros(1, size(ξ,2)); e[d] = h
    (f(ξ .+ e) - f(ξ .- e)) / 2h
end

@testset verbose = true "Core HighOrderMeshes.jl" begin

    @testset "Element Geometries" begin
        for D in 1:3, G in [Simplex, Block]
            G == Simplex && D == 1 && continue   # 1D meshes are always Block{1}
            geom = G{D}()

            @testset "$G{$D}" begin
                # 1. Basic properties
                @test dim(geom) == D
                @test nvertices(geom) > 0
                @test nfaces(geom) > 0
                @test nedges(geom) >= 0
                @test nnodes(geom, 1) == nvertices(geom)
                @test nnodes(geom, 3) == size(equispaced_nodes(geom, 3), 1)

                # 2. Consistency Checks
                fmap = facemap(geom)
                @test size(fmap) == (nvertices(subgeom(geom, D-1)), nfaces(geom))
                emap = edgemap(geom)
                @test size(emap) == (2, nedges(geom))

                # 3. Vertices and symmetries
                v = vertices(geom)
                @test size(v) == (nvertices(geom), D)
                @test all(0 .<= v .<= 1)
                syms = symmetries(geom)
                @test length(syms) == (G == Block ? 2^D * factorial(D) : factorial(D+1))
                vset = Set(Tuple.(eachrow(v)))
                for (A, b) in syms
                    @test Set(Tuple.(eachrow(v * A' .+ b'))) == vset
                end

                # 4. Inverse lookup and helper functions
                @test HOM.find_elgeom(D, nvertices(geom)) == geom
                @test subgeom(geom, D-1) isa ElementGeometry{D-1}
                @test length(HOM.midpoint(geom)) == D

                # 5. Show methods, check that it prints to a buffer
                s = sprint(show, geom)
                @test startswith(s, "ElementGeometry:")
                @test contains(s, "$(D)D")
            end
        end

        @test subgeom(Simplex{3}(), 1) == Block{1}()
        @test subgeom(Simplex{2}(), 0) == Block{0}()
        @test vertices(Block{2}())   == [0 0; 1 0; 0 1; 1 1]
        @test vertices(Simplex{3}()) == [0 0 0; 1 0 0; 0 1 0; 0 0 1]

        # Block faces: face 2d-1 is the side ξ_d = 0 and face 2d the side ξ_d = 1
        for D in 1:3
            v, fmap = vertices(Block{D}()), facemap(Block{D}())
            for d in 1:D, side in 0:1
                @test all(v[fmap[:, 2d-1+side], d] .== side)
            end
        end
        # Simplex faces: face i is opposite vertex i
        for D in 2:3
            fmap = facemap(Simplex{D}())
            for i in 1:D+1
                @test sort(fmap[:, i]) == setdiff(1:D+1, i)
            end
        end

        # Corner Cases / Error Handling
        @test_throws ErrorException HOM.find_elgeom(2, 100)
        @test_throws ErrorException HOM.plot_face_order(Simplex{4}())
        @test_throws ErrorException facemap(Simplex{4}())
        @test_throws ErrorException edgemap(Simplex{4}())
        @test_throws ErrorException facemap(Simplex{1}())
    end

    @testset "Generated geometry tables" begin
        # Hand-written tables before the extrusion-based generation, copied
        # verbatim as a regression check.
        @test facemap(Block{2}()) == [3 2 1 4; 1 4 2 3]
        @test facemap(Block{3}()) == [3 2 1 4 3 5; 1 4 2 3 4 6; 7 6 5 8 1 7; 5 8 6 7 2 8]

        @test edgemap(Block{2}()) == facemap(Block{2}())
        @test edgemap(Simplex{2}()) == facemap(Simplex{2}())
        @test edgemap(Block{3}()) == [1 3 5 7 1 2 5 6 1 2 3 4; 2 4 6 8 3 4 7 8 5 6 7 8]

        for D in 1:4
            v, fmap = vertices(Block{D}()), facemap(Block{D}())
            @test size(fmap) == (2^(D-1), 2D)
            for d in 1:D, side in 0:1
                @test all(v[fmap[:, 2d-1+side], d] .== side)
            end
        end

        # Orientation: det([n; a_1; ...; a_{D-1}]) > 0, where n is the outward
        # unit normal of the face and a_k = v_{k+1} - v_1 are its local axes
        # (v_1, ..., v_{2^{D-1}} being the face's vertex list).
        function block_face_sign(D, j)
            d, side = (j + 1) ÷ 2, (j - 1) % 2
            v, fmap = vertices(Block{D}()), facemap(Block{D}())
            vf = v[fmap[:, j], :]
            n  = zeros(D); n[d] = side == 0 ? -1 : 1
            axs = [vf[k+1, :] .- vf[1, :] for k in 1:D-1]
            sign(det(vcat(n', reduce(hcat, axs; init=zeros(D,0))')))
        end
        function simplex_face_sign(D, i)
            v, fmap = vertices(Simplex{D}()), facemap(Simplex{D}())
            vf  = v[fmap[:, i], :]
            n   = vec(sum(vf, dims=1))/size(vf,1) .- v[i, :]   # opposite vertex -> face centroid
            axs = [vf[k+1, :] .- vf[1, :] for k in 1:D-1]
            sign(det(vcat(n', reduce(hcat, axs; init=zeros(D,0))')))
        end
        for D in 2:3, j in 1:2D
            @test block_face_sign(D, j) > 0
        end
        for D in 2:3, i in 1:D+1
            @test simplex_face_sign(D, i) > 0
            @test sort(facemap(Simplex{D}())[:, i]) == setdiff(1:D+1, i)   # opposite vertex
        end

        @test nel(uniref(mshsquare(2), 1)) == 16
        let m = uniref(mshsquare(2), 1)
            v1, v2, v3, v4 = (m.x[m.el[k,:],:] for k in 1:4)
            area(a, b, c) = ((b[:,1].-a[:,1]).*(c[:,2].-a[:,2]) .- (c[:,1].-a[:,1]).*(b[:,2].-a[:,2]))/2
            @test all(area(v1,v2,v4) .+ area(v1,v4,v3) .> 0)
        end
        @test length(boundary_nodes(set_degree(uniref(mshsquare(2), 2), 2))) ==
              length(boundary_nodes(set_degree(mshsquare(8), 2)))
    end

    @testset "Polynomials" begin
        # Closed forms
        for x in (-0.7, 0.0, 0.3, 1.0)
            @test legendre(0, x) == 1
            @test legendre(1, x) ≈ x
            @test legendre(2, x) ≈ (3x^2 - 1)/2
            @test legendre(3, x) ≈ (5x^3 - 3x)/2
            @test dlegendre(3, x) ≈ (15x^2 - 3)/2
            @test jacobi(1, 1, 0, x) ≈ (3x + 1)/2
            @test legendre01(2, x) ≈ legendre(2, 2x-1)
            @test dlegendre01(2, x) ≈ 2dlegendre(2, 2x-1)
        end
        # Derivatives against finite differences
        for n in 0:5, (α, β) in ((0,0), (1,0), (3,0), (2,1))
            x, h = 0.37, 1e-6
            @test djacobi(n, α, β, x) ≈ (jacobi(n, α, β, x+h) - jacobi(n, α, β, x-h))/2h atol=1e-6
        end
        # Orthogonality of Jacobi polynomials with respect to (1-x)^α (1+x)^β
        x, w = gauss_legendre_quadrature(20)
        for (α, β) in ((0,0), (1,0), (3,0)), m in 0:4, n in 0:4
            val = sum(w .* (1 .- x).^α .* (1 .+ x).^β .* jacobi.(m, α, β, x) .* jacobi.(n, α, β, x))
            @test m == n ? val > 0 : abs(val) < 1e-12
        end
        # Exact arithmetic is preserved
        @test legendre(3, 1//2) == -7//16
        @test legendre(3, 1//2) isa Rational
        @test jacobi(2, 1, 0, 1//3) isa Rational

        # Node sets in the canonical order
        @test equispaced_nodes(Simplex{2}(), 2) == [0 0; 1//2 0; 1 0; 0 1//2; 1//2 1//2; 0 1]
        @test equispaced_nodes(Block{2}(), 1) == [0 0; 1 0; 0 1; 1 1]
        @test tensor_nodes(Block{2}(), [0.0, 0.5, 1.0]) == Float64.(equispaced_nodes(Block{2}(), 2))
        @test size(equispaced_nodes(Simplex{3}(), 4)) == (35, 3)
        @test size(equispaced_nodes(Block{0}(), 3)) == (1, 0)
        @test multiindices(Simplex{2}(), 2) == [(0,0), (1,0), (2,0), (0,1), (1,1), (0,2)]

        # Polynomial bases: sizes and gradients
        for eg in (Simplex{2}(), Simplex{3}(), Block{1}(), Block{2}(), Block{3}()), p in (1, 3, 5)
            D = dim(eg)
            ξ = Float64.(equispaced_nodes(eg, p)) .* 0.8 .+ 0.05
            V  = polybasis(eg, ξ, p)
            dV = dpolybasis(eg, ξ, p)
            @test size(V)  == (size(ξ,1), nnodes(eg, p))
            @test size(dV) == (size(ξ,1), nnodes(eg, p), D)
            for d in 1:D
                fd = fdgrad(s -> polybasis(eg, s, p), ξ, d)
                @test maximum(abs.(dV[:,:,d] - fd)) < 1e-6 * max(1, maximum(abs.(fd)))
            end
        end
        # The simplex basis is orthogonal and finite at the vertices
        for (eg, p) in ((Simplex{2}(), 6), (Simplex{3}(), 4))
            ξ, w = quadrature(eg, 2p)
            V = polybasis(eg, ξ, p)
            M = V' * (w .* V)
            @test maximum(abs.(M - Diagonal(M))) < 1e-10 * maximum(abs.(M))
            v = Float64.(vertices(eg))
            @test all(isfinite, polybasis(eg, v, 4))
            @test all(isfinite, dpolybasis(eg, v, 4))
        end
        @test polybasis(Block{0}(), zeros(3, 0), 2) == ones(3, 1)

        # Quadrature: exact for all monomials of degree <= p, for both geometries
        for D in 1:3, p in 1:12
            ξ, w = quadrature(Block{D}(), p)
            @test length(w) == cld(p+1, 2)^D
            for α in multiindices(Block{D}(), p)
                sum(α) <= p || continue
                @test w' * monomial(ξ, α) ≈ prod(1 ./ (α .+ 1))
            end
        end
        for eg in (Simplex{2}(), Simplex{3}()), p in 1:HOM.max_quadrature_degree(eg)
            D = dim(eg)
            ξ, w = quadrature(eg, p)
            for α in multiindices(eg, p)
                exact = Float64(prod(factorial.(big.(α))) // factorial(big(sum(α) + D)))
                @test w' * monomial(ξ, α) ≈ exact atol=1e-11   # tabulated rules are accurate to ~1e-12
            end
        end
        @test eltype(quadrature(Block{2}(), 3, Float32)[2]) == Float32
        @test eltype(quadrature(Simplex{2}(), 3, Float32)[1]) == Float32
        @test_throws ErrorException quadrature(Simplex{2}(), 1000)
        @test Base.return_types(quadrature, (Simplex{2}, Int))[1] == Tuple{Matrix{Float64},Vector{Float64}}

        # 1D rules
        x,w = gauss_legendre_quadrature(3)
        @test w' * x.^4 ≈ 2/5
        x,w = gauss_lobatto_quadrature(4)
        @test w' * x.^4 ≈ 2/5
        x,w = gauss_legendre01_quadrature(3)
        @test w' * x.^4 ≈ 1/5
        x,w = gauss_lobatto01_quadrature(4)
        @test w' * x.^4 ≈ 1/5
        @test gauss_lobatto01_nodes(5) ≈ 1 .- reverse(gauss_lobatto01_nodes(5))
        @test eltype(gauss_legendre_nodes(4, Float32)) == Float32
    end

    @testset "Reference element" begin
        for eg in (Simplex{2}(), Simplex{3}(), Block{1}(), Block{2}(), Block{3}()), p in (1, 2, 4)
            fe = FiniteElement(eg, p)
            D  = dim(eg)
            @test porder(fe) == p && nnodes(fe) == nnodes(eg, p)
            @test elgeom(fe) == eg && dim(fe) == D

            # Kronecker delta at the nodes, partition of unity elsewhere
            @test shapefcns(fe, ref_nodes(fe)) ≈ I atol=1e-10
            ξ  = Float64.(equispaced_nodes(eg, p+2)) .* 0.8 .+ 0.05
            N  = shapefcns(fe, ξ)
            dN = dshapefcns(fe, ξ)
            @test all(isapprox.(sum(N, dims=2), 1, atol=1e-10))
            @test all(abs.(sum(dN, dims=2)) .< 1e-8)
            for d in 1:D
                fd = fdgrad(s -> shapefcns(fe, s), ξ, d)
                @test maximum(abs.(dN[:,:,d] - fd)) < 1e-6 * max(1, maximum(abs.(fd)))
            end

            # Corner nodes, sub-elements and restrictions
            @test ref_nodes(fe)[corner_nodes(fe), :] == vertices(eg)
            for d in 0:D
                sub = subelement(fe, d)
                @test elgeom(sub) == subgeom(eg, d) && porder(sub) == p
                @test ref_nodes(sub) == ref_nodes(fe, d)
                @test ref_nodes(sub) == Float64.(equispaced_nodes(subgeom(eg, d), p))
            end

            # The explicit constructor infers the degree and accepts the node set
            @test porder(FiniteElement(eg, ref_nodes(fe))) == p

            # Interpolation
            u = randn(nnodes(fe), 3)
            @test interpolate(shapefcns(fe, ref_nodes(fe)), u) ≈ u
            @test interpolate(N, u[:,1]) ≈ N * u[:,1]
            @test size(interpolate(dN, u)) == (size(ξ,1), 3, D)
            @test size(interpolate(N, randn(nnodes(fe), 3, 2))) == (size(ξ,1), 3, 2)
            @test interpolate(dN, u)[:,:,D] ≈ dN[:,:,D] * u
        end

        # Gauss-Lobatto tensor elements
        for D in 1:3, p in (2, 5)
            fe = FiniteElement(Block{D}(), gauss_lobatto01_nodes(p+1))
            @test porder(fe) == p
            @test shapefcns(fe, ref_nodes(fe)) ≈ I atol=1e-10
            @test ref_nodes(fe, 1)[:, 1] ≈ gauss_lobatto01_nodes(p+1)
        end

        # Non-conforming node sets are rejected
        @test_throws ArgumentError FiniteElement(Block{2}(), [0.0, 0.3, 1.0])      # non-symmetric line
        @test_throws ArgumentError FiniteElement(Block{1}(), [0.0, 0.5, 0.9])      # missing vertex
        bad = Float64.(equispaced_nodes(Simplex{2}(), 3)); bad[2, 1] += 0.01        # edge node moved
        @test_throws ArgumentError FiniteElement(Simplex{2}(), bad)
        @test_throws ArgumentError FiniteElement(Simplex{2}(), zeros(7, 2))        # no matching degree
        @test_throws ArgumentError FiniteElement(Simplex{2}(), zeros(6, 3))        # wrong dimension
        @test FiniteElement(Simplex{2}(), bad; check=false) isa FiniteElement
        @test_throws ArgumentError check_conforming(Block{2}(), tensor_nodes(Block{2}(), [0.0, 0.3, 1.0]))
        @test isnothing(check_conforming(Block{3}(), tensor_nodes(Block{3}(), gauss_lobatto01_nodes(4))))

        # Linear shape functions on a geometry
        @test shapefcns(Simplex{2}(), [0.25 0.25]) ≈ [0.5 0.25 0.25]
        @test shapefcns(Block{2}(), [0.5 0.5]) ≈ [0.25 0.25 0.25 0.25]
        @test dshapefcns(Block{1}(), [0.3]) ≈ reshape([-1.0 1.0], 1, 2, 1)

        # Exact arithmetic
        fe = FiniteElement(Block{2}(), 2, Rational{Int})
        @test eltype(ref_nodes(fe)) == Rational{Int}
        @test shapefcns(fe, ref_nodes(fe)) == I

        @test startswith(sprint(show, FiniteElement(Simplex{2}(), 3)), "FiniteElement: p=3 triangle")

        # Type stability of the core
        fe = FiniteElement(Block{2}(), 3); ξ = ref_nodes(fe)
        @test Base.return_types(shapefcns, (typeof(fe), typeof(ξ)))[1] == Matrix{Float64}
        @test Base.return_types(dshapefcns, (typeof(fe), typeof(ξ)))[1] == Array{Float64,3}
        @test Base.return_types(legendre, (Int, Float64))[1] == Float64
        @test Base.return_types(quadrature, (Block{2}, Int))[1] == Tuple{Matrix{Float64},Vector{Float64}}
        m = mshsquare(2)
        @test Base.return_types(set_degree, (typeof(m), Int))[1] == typeof(m)
    end

    @testset "Basic ex1mesh properties" begin
        for eg in (Block{2}(), Simplex{2}()), nref in 1:3
            m = ex1mesh(nref=nref, eg=eg)
            @test porder(m) == 3
            @test elgeom(m) == eg
            @test m isa HighOrderMesh{2,typeof(eg),Float64}
        end
    end

    @testset "Basic mshsquare properties" begin
        for xtype in (Float64, Float32, Rational{Int})
            m = mshsquare(5, 3; T=xtype)
            @test size(m.el) == (4,15)
            @test eltype(m.x) == xtype
        end
        @test contains(sprint(show, mshsquare(2)), "9 nodes")
    end

    @testset "Mesh generators" begin
        # Simplex hypercubes: conforming, positively oriented, unit volume, faces tagged by plane
        for (eg, dims, p) in ((Simplex{2}(), (3, 2), 3), (Simplex{3}(), (2, 2, 1), 2))
            D = length(dims)
            m = mshhypercube(dims; eg, p)
            @test elgeom(m) == eg
            @test nel(m) == prod(dims) * factorial(D)
            @test nnodes(m) == prod(p .* dims .+ 1)
            vol, q = volume_and_quality(m)
            @test vol ≈ 1 && q ≈ 1
            for d in 1:D, (tag, c) in ((2d-1, 0), (2d, 1))
                X = boundary_face_nodes(m, tag)
                @test length(X) == prod(dims[[1:d-1; d+1:D]]) * factorial(D-1)
                @test all(all(abs.(X[:,d] .- c) .< 1e-14) for X in X)
            end
        end
        @test mshsquare(8, 4; eg=Simplex{2}(), p=3) isa HighOrderMesh{2,Simplex{2}}
        @test nel(mshcube(4; eg=Simplex{3}())) == 384
        @test_throws ArgumentError mshsquare(2; eg=Simplex{3}())
        for m in (mshsquare(3; eg=Simplex{2}(), periodic_dirs=(1,2)),
                  mshcube(2; eg=Simplex{3}(), periodic_dirs=(1,2,3)))
            @test !any(isboundary, m.nb)
            @test all(m.nb[nb.face,nb.el].el == iel for iel in axes(m.nb,2) for nb in m.nb[:,iel])
        end

        # Triangle disks: nodes on the arcs exactly on the circle and at equal angles
        for (shape, nsec, β) in ((:full, 6, 1/3), (:half, 3, 1/3), (:quarter, 2, 1/4)), n in 1:3, p in 1:3
            m = mshcircle(n; eg=Simplex{2}(), p, shape)
            @test nel(m) == nsec * n^2
            arc = maximum(bndtag.(filter(isboundary, m.nb)))
            X = reduce(vcat, boundary_face_nodes(m, arc))
            @test all(abs.(hypot.(X[:,1], X[:,2]) .- 1) .< 1e-14)
            θ = sort(mod.(atan.(X[:,2], X[:,1]) / π, 2))
            θ = θ[[true; diff(θ) .> 1e-9]]
            @test θ ≈ (0:nsec*n*p - (shape == :full)) .* (β/(n*p)) atol=1e-12
            shape == :quarter && @test all(all(X[:,2] .== 0) for X in boundary_face_nodes(m, 1)) &&
                                       all(all(abs.(X[:,1]) .< 1e-15) for X in boundary_face_nodes(m, 2))
        end
        for p in 1:4
            errs = [ abs(volume_and_quality(mshcircle(n; eg=Simplex{2}(), p))[1] - π) for n in (2, 4) ]
            @test errs[1] / errs[2] > 2^(p+1) * 0.9    # area converges at least like h^(p+1)
            @test volume_and_quality(mshcircle(4; eg=Simplex{2}(), p))[2] > 0.8
        end
        @test_throws ArgumentError mshcircle(2; shape=:third)

        # Gmsh samples
        @test_throws ArgumentError gmsh_sample(:torus)
        @test_throws ArgumentError gmsh_sample(:circle; eg=Simplex{3}())
        if Sys.which("gmsh") !== nothing
            for eg in (Simplex{2}(), Block{2}())
                m = gmsh_sample(:square; h=0.3, eg)
                @test elgeom(m) == eg
                @test volume_and_quality(m)[1] ≈ 1
                for d in 1:2, (tag, c) in ((2d-1, 0), (2d, 1))
                    @test all(all(abs.(X[:,d] .- c) .< 1e-14) for X in boundary_face_nodes(m, tag))
                end
                m = gmsh_sample(:circle; h=0.3, p=3, eg)
                @test elgeom(m) == eg
                @test all(abs.(norm.(eachrow(m.x[boundary_nodes(m),:])) .- 1) .< 1e-12)
                @test volume_and_quality(m)[1] ≈ π rtol=1e-4
                @test volume_and_quality(m)[2] > 0.2
            end
            for h in (0.2, 0.1)    # no nearly flat quad corners at the boundary
                @test volume_and_quality(gmsh_sample(:circle; h, p=3, eg=Block{2}()))[2] > 0.2
            end
            for eg in (Simplex{3}(), Block{3}())
                m = gmsh_sample(:cube; h=0.5, eg)
                @test elgeom(m) == eg
                @test volume_and_quality(m)[1] ≈ 1
                @test sort(unique(bndtag.(filter(isboundary, m.nb)))) == 1:6
                @test all(all(abs.(X[:,3] .- 1) .< 1e-14) for X in boundary_face_nodes(m, 6))
            end
            m = gmsh_sample(:sphere; h=0.5, p=2)
            @test volume_and_quality(m)[1] ≈ 4π/3 rtol=1e-2
        end
    end

    @testset "Neighbor" begin
        @test Neighbor(2, 1, 3) == Neighbor(2, 1, 3)
        @test isboundary(Neighbor(-2, 0, 0)) && bndtag(Neighbor(-2, 0, 0)) == 2
        @test !isboundary(Neighbor(5, 1, 1))
        @test sizeof(Neighbor) == 8
        m = mshsquare(2)
        @test count(isboundary, m.nb) == 8
    end

    @testset "Degree and node-set changes" begin
        # Changing the degree up and down reproduces the same mesh
        for m in (mshsquare(2), ex1mesh(eg=Simplex{2}(), nref=1))
            m2  = set_degree(m, 2)
            m42 = set_degree(set_degree(m, 4), 2)
            @test m42.x ≈ m2.x atol=1e-12
            @test m42.el == m2.el
            @test m42.nb == m2.nb
            @test m2.nb !== m.nb          # no sharing between meshes
        end
        let m = ex1mesh(eg=Simplex{2}(), nref=1)
            m3  = set_degree(m, 3)
            m43 = set_degree(set_degree(m, 4), 3)
            @test m43.x ≈ m3.x atol=1e-12
            @test m43.el == m3.el
            @test m43.nb == m3.nb
        end
        # Lobatto nodes on blocks, and an explicit node matrix
        m  = set_degree(mshcube(2, 1, 1), 3)
        ml = set_lobatto_nodes(m)
        @test porder(ml) == 3 && nnodes(ml) == nnodes(m)
        @test ref_nodes(ml.fe, 1)[:,1] ≈ gauss_lobatto01_nodes(4)
        me = set_ref_nodes(m, ref_nodes(ml.fe))
        @test me.x ≈ ml.x
        @test_throws ArgumentError set_ref_nodes(m, [0.0, 0.2, 0.5, 1.0])
        # Shared faces carry the same node sets in 3D
        for mm in (set_degree(mshcube(2,2,2), 3), set_lobatto_nodes(set_degree(mshcube(2,2,2), 4)))
            f2n = mkface2nodes(mm)
            for iel in 1:nel(mm), j in 1:nfaces(elgeom(mm))
                jel, k = mm.nb[j,iel].el, mm.nb[j,iel].face
                jel > 0 || continue
                @test Set(mm.el[f2n[:,j],iel]) == Set(mm.el[f2n[:,k],jel])
            end
        end
    end

    @testset "IO - .hom format" begin
        rootdir = pkgdir(HighOrderMeshes)
        meshes = HighOrderMesh[ex1mesh(nref=nref, eg=eg) for eg in (Block{2}(), Simplex{2}()) for nref in 0:1]
        push!(meshes, set_lobatto_nodes(set_degree(mshcube(2,1,1), 3)))
        push!(meshes, mshsquare(3; T=Float32))
        for m in meshes
            mktempdir() do tmpdir
                tmp_path = joinpath(tmpdir, "test_mesh.hom")
                savemesh(tmp_path, m)
                m2 = loadmesh(tmp_path)
                @test typeof(m2) == typeof(m)
                @test m.x == m2.x && m.el == m2.el && m.nb == m2.nb
                @test ref_nodes(m.fe) == ref_nodes(m2.fe)
                @test porder(m2) == porder(m)
            end
        end
        m = loadmesh(joinpath(rootdir, "examples", "meshes", "naca_mesh_1.hom"))
        @test m isa HighOrderMesh
        mktempdir() do tmpdir
            tmp_path = joinpath(tmpdir, "bad.hom")
            write(tmp_path, "NOTAMESH")
            @test_throws ErrorException loadmesh(tmp_path)
        end
    end

    @testset "gmsh file import" begin
        rootdir = pkgdir(HighOrderMeshes)
        for (filename,eg) in (("circle_tris.msh", Simplex{2}()),
                              ("square_tris.msh", Simplex{2}()),
                              ("circle_quads.msh", Block{2}()),
                              ("square_quads.msh", Block{2}()))
            fullname = joinpath(rootdir, "examples", "gmsh", filename)
            m = gmsh2msh(fullname; verbose=false)
            @test elgeom(m) == eg
        end

        circle = joinpath(rootdir, "examples", "gmsh", "circle_tris.msh")
        @test gmsh_physical_names(circle) == Dict(1 => "Circle")
        let fname = tempname() * ".msh"          # physical names may contain spaces
            write(fname, replace(read(circle, String), "\"Circle\"" => "\"Unit circle\""))
            @test gmsh_physical_names(fname) == Dict(1 => "Unit circle")
            rm(fname)
        end
        m_verbose = gmsh2msh(circle)
        m_quiet   = gmsh2msh(circle; verbose=false)
        @test m_verbose.x == m_quiet.x && m_verbose.el == m_quiet.el && m_verbose.nb == m_quiet.nb
    end

    @testset "Gmsh affine round trip" begin
        if Sys.which("gmsh") !== nothing
            # A straight-sided mesh is the affine (tris, tets) or bilinear
            # (transfinite quads, hexes) image of the reference element, so
            # the high-order nodes must equal the p=1 interpolation of the
            # corner nodes.
            affine_file = joinpath(pkgdir(HighOrderMeshes), "docs", "dev", "gmsh_affine_test.jl")
            if isfile(affine_file)
                include(affine_file)
                geos = (geo_square, geo_square_quads, geo_box, geo_box_hexes)
            else
                fallback_square = """
                    Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
                    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
                    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
                    """
                fallback_square_quads = fallback_square *
                    "Transfinite Curve{1,2,3,4}=4; Transfinite Surface{1}; Recombine Surface{1};\n"
                fallback_box = """
                    SetFactory("OpenCASCADE");
                    Box(1)={0,0,0,1,1,1};
                    """
                fallback_box_hexes = """
                    Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
                    Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
                    Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
                    Transfinite Curve{1,2,3,4}=3; Transfinite Surface{1}; Recombine Surface{1};
                    Extrude {0,0,1} { Surface{1}; Layers{2}; Recombine; };
                    """
                geos = (fallback_square, fallback_square_quads, fallback_box, fallback_box_hexes)
            end

            for geo in geos, p in 1:5
                m       = gmshstr2msh(geo; p, cmdadd="-v 0")
                fe      = m.fe
                xdg     = dg_nodes(m)
                corners = xdg[corner_nodes(fe), :, :]
                N1      = shapefcns(elgeom(m), ref_nodes(fe))
                pred    = interpolate(N1, corners)
                @test maximum(abs.(pred - xdg)) < 1e-10
            end
        end
    end

    @testset "FEM assembly" begin
        # L2 error of a DG field U against the function uex, with the quadrature of d
        l2err(d, U, uex) = sqrt(sum(d.met.wJ .* (d.ref.ϕ * U .- uex.(d.met.x)).^2))
        rate(errs) = log2(errs[1] / errs[2])

        @testset "CG Poisson" begin
            # Solve -∇²u = 1 with zero Dirichlet boundary conditions on the unit circle
            for n = 1:4, porder = 1:4
                m = mshcircle(n, p=porder)
                d = CGData(m)
                A = assemble_matrix(m.el, elmats(elmat_laplace!, d))
                f = assemble_vector(m.el, elvec_source(d, xy -> 1.0))
                A, f = strong_dirichlet(A, f, boundary_nodes(m))
                u = A \ f
                uexact = (1 .- sum(m.x.^2, dims=2)) / 4
                error = maximum(abs.(u[:] - uexact[:]))
                # Assume O(h^{p+1}) convergence, with fitted constant (upper bound)
                error_bound = (0.2 / n) ^ (porder + 1)
                @test error < error_bound
            end
        end

        @testset "CG mass, residual and Dirichlet data" begin
            m = mshcircle(2, p=3)
            d = CGData(m)
            M = assemble_matrix(m.el, elmats(elmat_mass!, d))
            @test sum(M) ≈ π rtol=1e-4                       # area of the (p=3) unit disk
            A = assemble_matrix(m.el, elmats(elmat_laplace!, d))
            @test norm(A - A') < 1e-12 * norm(A)
            u0 = randn(nnodes(m))
            R  = laplace_residual!(zeros(size(m.el)), zeros(2npoints(d.ref), nel(m)), d, u0[m.el])
            @test assemble_vector(m.el, R) ≈ A * u0            # matrix-free residual
            bnd = boundary_nodes(m)
            A, f = strong_dirichlet(A, zeros(nnodes(m)), bnd, m.x[bnd, 1])
            @test A \ f ≈ m.x[:, 1] atol=1e-10                # u = x is harmonic and in the space
        end

        @testset "DG interior penalty" begin
            uex(x) = sin(pi*x[1]) * sin(pi*x[2]);  fex(x) = 2pi^2 * uex(x)
            for p in 1:3
                errs = Float64[]
                for n in (4, 8)
                    m = mshsquare(n, p=p)
                    d = DGData(m)
                    A, b = dg_laplace(d, Dict(k => uex for k in 1:4))
                    @test norm(A - A') < 1e-12 * norm(A)
                    U = reshape(A \ vec(elvec_source(d, fex) + b), size(m.el))
                    push!(errs, l2err(d, U, uex))
                end
                @test rate(errs) > p + 0.5
            end
            # curved triangles, Dirichlet on boundaries 1 and 3, Neumann on 2 and 4
            uex2(x) = cos(pi*x[1]) * cos(pi*x[2]);  fex2(x) = 2pi^2 * uex2(x)
            errs = Float64[]
            for nref in (2, 3)
                m = ex1mesh(eg=Simplex{2}(), nref=nref)
                d = DGData(m)
                A, b = dg_laplace(d, Dict(1 => uex2, 3 => uex2))
                U = reshape(A \ vec(elvec_source(d, fex2) + b), size(m.el))
                push!(errs, l2err(d, U, uex2))
            end
            @test rate(errs) > 3.4
        end

        @testset "DG convection and convection-diffusion" begin
            vel(x) = (1.0, 2x[1])                                  # divergence free
            m = mshsquare(6, p=3)
            d = DGData(m)
            A, b = dg_convection(d, vel, Dict(k => x -> 1.0 for k in 1:4))
            @test maximum(abs, A \ vec(b) .- 1) < 1e-10           # constant state preserved
            uex(x) = sin(pi*x[1]) * sin(pi*x[2])
            fex(x) = pi*cos(pi*x[1])*sin(pi*x[2]) + 2x[1]*pi*sin(pi*x[1])*cos(pi*x[2])
            for p in (1, 3)
                errs = Float64[]
                for n in (4, 8)
                    m = mshsquare(n, p=p)
                    d = DGData(m)
                    A, b = dg_convection(d, vel, Dict(1 => uex, 3 => uex))   # inflow boundaries
                    U = reshape(A \ vec(elvec_source(d, fex) + b), size(m.el))
                    push!(errs, l2err(d, U, uex))
                end
                @test rate(errs) > p + 0.5
            end
            ε = 0.1
            uex3(x) = sin(pi*x[1]/2) * sin(pi*x[2]/2)              # ∂u/∂n = 0 on the outflow boundaries
            fcd(x) = ε*(pi^2/2)*uex3(x) + (pi/2)*cos(pi*x[1]/2)*sin(pi*x[2]/2) + 2x[1]*(pi/2)*sin(pi*x[1]/2)*cos(pi*x[2]/2)
            errs = Float64[]
            for n in (4, 8)
                m = mshsquare(n, p=3)
                d = DGData(m)
                bc = Dict(1 => uex3, 3 => uex3)
                C, bC = dg_convection(d, vel, bc)
                K, bK = dg_laplace(d, bc)
                U = reshape((C + ε*K) \ vec(elvec_source(d, fcd) + bC + ε*bK), size(m.el))
                push!(errs, l2err(d, U, uex3))
            end
            @test rate(errs) > 3.5
        end
    end

    @testset "Converters / export" begin
        ## VTK
        for eg in (Simplex{2}(), Block{2}())
            m = ex1mesh(eg=eg)

            mktempdir() do tmpdir
                u = ex1solution(m)
                filename = joinpath(tmpdir, "ex1.vtk")
                vtkwrite(filename, m, u)

                @test isfile(filename)
                @test filesize(filename) > 0
                header = open(readline, filename)
                @test startswith(header, "# vtk")
            end

            mktempdir() do tmpdir
                u = hcat(m.x[:,1].^2, m.x[:,2], -m.x[:,1])
                filename = joinpath(tmpdir, "ex1.vtk")
                vtkwrite(filename, m, u, umap=[1, 2:3])

                @test isfile(filename)
                @test filesize(filename) > 0
                header = open(readline, filename)
                @test startswith(header, "# vtk")
            end
        end

        ## 3DG export
        for eg in (Simplex{2}(), Block{2}())
            m = ex1mesh(eg=eg)
            if eg == Block{2}()
                m = set_lobatto_nodes(m)
            end

            flds = mshto3dg(m)
            @test size(flds.p1)[[1,3]] == size(m.el)

            m2 = mshfrom3dg(flds)
            @test elgeom(m2) == elgeom(m)
            @test porder(m2) == porder(m)
            @test dg_nodes(m2) ≈ dg_nodes(m)
            @test m2.nb == m.nb

            flds2 = mshto3dg(m2)
            @test flds2.p1 ≈ flds.p1
            @test flds2.t2t == flds.t2t
            @test flds2.t2n == flds.t2n
        end

        m = mshcube(2, 1, 1; p=1)
        flds = mshto3dg(m)
        t2n = copy(flds.t2n)
        i = findfirst(>=(0), flds.t2t)
        perm = 5
        t2n[i] = (t2n[i] & 0x0f) + 2^7 + 2^4 * (perm - 1)
        m2 = mshfrom3dg(merge(flds, (; t2n)))
        @test m2.nb[i].face == (flds.t2n[i] & 0x0f) + 1
        @test m2.nb[i].perm == perm
    end

    @testset "Mesh utilities" begin
        msh = mshsquare(5)
        set_bnd_periodic!(msh, (1,2), 1)    # Periodic left/right (x-direction)
        set_bnd_periodic!(msh, (3,4), 2)    # Periodic bottom/top (y-direction)
        @test minimum(getfield.(msh.nb, :el)[:]) == 1  # No actual boundaries
        @test all(getfield.(msh.nb, :perm)[:] .== 1)

        x = [0.0 0.0 0.0
             1.0 0.0 0.0
             0.0 1.0 0.0
             0.0 0.0 1.0
             1.0 1.0 1.0]
        el = [1 5
              2 3
              3 4
              4 2]
        nb = HighOrderMeshes.el2nb(el, Simplex{3}())
        @test nb[1,1] == Neighbor(2, 1, 3)
        @test nb[1,2] == Neighbor(1, 1, 2)

        msh = mshcube(2, 1, 1)
        oldnb = copy(msh.nb)
        set_bnd_periodic!(msh, (1,2), 1)
        periodic_faces = findall(i -> oldnb[i].el < 0 && msh.nb[i].el > 0, eachindex(msh.nb))
        @test !isempty(periodic_faces)
        @test all(i -> msh.nb[i].perm >= 1, periodic_faces)

        msh = ex1mesh(nref=1, eg=Block{2}())
        align_with_ldgswitch!(msh)
        sw = mkldgswitch(msh)
        @test all(sw[1,:] .== 1 .&& sw[3,:] .== 1)

        @test nel(uniref(mshsquare(2), 2)) == 64
        @test length(boundary_nodes(set_degree(mshsquare(2), 3))) == 24
        @test length(boundary_nodes(set_degree(mshsquare(2), 3), 1)) == 7
    end

    @testset "Boundary distance" begin
        # Straight faces: exact. mshsquare numbers x=0 as 1 and y=0 as 3.
        m = mshsquare(3; p=3)
        x = dg_nodes(m)
        X, Y = x[:,:,1], x[:,:,2]
        @test size(boundary_distance(m)) == size(X)
        @test boundary_distance(m, 1) ≈ X atol=1e-14
        @test boundary_distance(m, (1,3)) ≈ min.(X, Y) atol=1e-14
        @test boundary_distance(m) ≈ min.(X, Y, 1 .- X, 1 .- Y) atol=1e-14
        @test_throws ErrorException boundary_distance(m, 7)

        cubedist(x) = min.(x[:,:,1], x[:,:,2], x[:,:,3], 1 .- x[:,:,1], 1 .- x[:,:,2], 1 .- x[:,:,3])
        m = mshcube(2, 2, 2; p=2)
        @test boundary_distance(m; nsub=4) ≈ cubedist(dg_nodes(m)) atol=1e-14

        # Unit cube split into 6 tetrahedra (triangular faces)
        vix(v) = 1 + v[1] + 2v[2] + 4v[3]
        el = Int[]
        for σ in ((1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1))
            v = [0, 0, 0]
            push!(el, vix(v))
            for k in σ
                v[k] += 1
                push!(el, vix(v))
            end
        end
        xv = Float64[ v[c] for v in [ (i,j,k) for k in 0:1 for j in 0:1 for i in 0:1 ], c in 1:3 ]
        m = set_degree(HighOrderMesh(xv, reshape(el, 4, 6); bndexpr=x -> [x'; 1 .- x'][:]), 3)
        x = dg_nodes(m)
        @test boundary_distance(m; nsub=4) ≈ cubedist(x) atol=1e-14
        @test boundary_distance(m, 1; nsub=4) ≈ x[:,:,1] atol=1e-14

        # Curved faces: close to the exact distance, and O(h^2) convergence in nsub
        m = mshcircle(2; p=4)
        x = dg_nodes(m)
        r = sqrt.(x[:,:,1].^2 .+ x[:,:,2].^2)
        d8, d32, dref = (boundary_distance(m; nsub=n) for n in (8, 32, 256))
        @test maximum(abs.(d32 .- (1 .- r))) < 1e-4
        @test maximum(abs.(d8 .- dref)) / maximum(abs.(d32 .- dref)) > 10

        if Sys.which("gmsh") !== nothing
            m = gmsh_sample(:sphere; h=0.3, p=2)
            x = dg_nodes(m)
            r = sqrt.(sum(x.^2, dims=3)[:,:,1])
            d4, d16, dref = (boundary_distance(m; nsub=n) for n in (4, 16, 64))
            @test maximum(abs.(d16 .- (1 .- r))) < 1e-3   # p=2 geometry error ~2.5e-4
            @test maximum(abs.(d4 .- dref)) / maximum(abs.(d16 .- dref)) > 10
        end
    end

    @testset "Refinement" begin
        # Area by quadrature of the Jacobian determinant; errors on inverted elements.
        function mesh_area(m)
            ξ, w = quadrature(elgeom(m), 3*porder(m) + 2)
            J = interpolate(dshapefcns(m.fe, ξ), dg_nodes(m))
            dets = J[:,:,1,1] .* J[:,:,2,2] .- J[:,:,1,2] .* J[:,:,2,1]
            @assert minimum(dets) > 0
            sum(w .* dets)
        end
        # Node count of a conforming mesh: vertices, edge nodes, interior nodes.
        function conforming_nnodes(m)
            p   = porder(m)
            nf  = size(m.nb, 1)
            nv  = length(unique(m.el[corner_nodes(m.fe),:]))
            ne  = (count(isboundary, m.nb) + nf*nel(m)) ÷ 2
            nint = size(m.el, 1) - nvertices(elgeom(m)) - nf*(p - 1)
            nv + (p - 1)*ne + nint*nel(m)
        end
        # Every boundary face has its nodes on the curve of its tag.
        function tags_consistent(m, bndexpr; tol=1e-12)
            f2n = mkface2nodes(m)
            all(isboundary(m.nb[j,iel]) ?
                all(abs(bndexpr(m.x[k,:])[bndtag(m.nb[j,iel])]) < tol for k in m.el[f2n[:,j],iel]) : true
                for j in axes(m.nb,1), iel in axes(m.nb,2))
        end
        # Unit square with 4 boundary tags and a curved bottom edge, as quads or triangles.
        # The map is a cubic, so p=3 elements and their children represent it exactly.
        sqbnd(x) = [x[1], 1 - x[1], x[2] - 0.8*x[1]*(1 - x[1]), 1 - x[2]]
        function curved_square(n, eg; p=3)
            m = mshsquare(n)
            if eg isa Simplex
                el = hcat(m.el[[1,2,4],:], m.el[[1,4,3],:])
                m  = HighOrderMesh(m.x, el; bndexpr=x -> [x[1], 1 - x[1], x[2], 1 - x[2]])
            end
            m = set_degree(m, p)
            m.x[:,2] .+= 0.8*m.x[:,1].*(1 .- m.x[:,1]).*(1 .- m.x[:,2])
            m
        end

        # Uniform refinement of curved meshes reproduces the geometry exactly.
        for eg in (Block{2}(), Simplex{2}())
            for m in (ex1mesh(nref=2, eg=eg), curved_square(2, eg))
                m2 = uniref(m, 2)
                @test nel(m2) == 16*nel(m)
                @test mesh_area(m2) ≈ mesh_area(m) rtol=1e-13
                @test nnodes(m2) == conforming_nnodes(m2)
                @test count(isboundary, m2.nb) == 4*count(isboundary, m.nb)
            end
            m = curved_square(2, eg)
            @test tags_consistent(uniref(m, 2), sqbnd)
            @test sort(unique(bndtag.(filter(isboundary, uniref(m).nb)))) == 1:4
        end
        @test sortslices(uniref(mshsquare(2)).x, dims=1) == sortslices(mshsquare(4).x, dims=1)
        @test nnodes(uniref(ex1mesh(nref=1), 0)) == nnodes(ex1mesh(nref=1))

        # Single-element patterns: marked edges => number of children (after closure).
        mq = curved_square(1, Block{2}())
        a0 = mesh_area(mq)
        # The area depends only on the boundary curves, so it is exact for every pattern.
        for (marks, nchildren) in ((Bool[0,0,1,1], 2), (Bool[1,1,0,0], 2),
                                   (Bool[1,0,1,0], 3), (Bool[0,1,0,1], 3),
                                   (Bool[1,0,0,1], 3), (Bool[0,1,1,0], 3),
                                   (Bool[1,0,0,0], 3), (Bool[1,1,1,0], 4),
                                   (Bool[1,1,1,1], 4), (Bool[0,0,0,0], 1))
            m1, parent = refine_with_parents(mq, reshape(marks, 4, 1))
            @test nel(m1) == nchildren
            @test parent == ones(Int, nchildren)
            @test mesh_area(m1) ≈ a0 rtol=1e-13
            @test nnodes(m1) == conforming_nnodes(m1)
            @test tags_consistent(m1, sqbnd)
        end
        mt = HighOrderMesh(Float64[0 0; 1 0; 0 1], reshape([1, 2, 3], 3, 1);
                           bndexpr=x -> [x[1] + x[2] - 1, x[1], x[2]])
        mt = set_degree(mt, 3)
        mt.x[:,1] .+= 0.1*mt.x[:,2].*(1 .- mt.x[:,2])              # curves faces 1 and 2
        tribnd(x) = [x[1] - 0.1*x[2]*(1 - x[2]) + x[2] - 1, x[1] - 0.1*x[2]*(1 - x[2]), x[2]]
        for (marks, nchildren) in ((Bool[1,0,0], 2), (Bool[0,1,0], 2), (Bool[0,0,1], 2),
                                   (Bool[1,1,0], 4), (Bool[1,1,1], 4))
            m1 = refine(mt, reshape(marks, 3, 1))
            @test nel(m1) == nchildren
            @test mesh_area(m1) ≈ mesh_area(mt) rtol=1e-13
            @test nnodes(m1) == conforming_nnodes(m1)
            @test tags_consistent(m1, tribnd)
        end

        # Marks propagate through the closure and stay conforming.
        for m in (set_degree(mshsquare(4), 2), ex1mesh(nref=2, eg=Simplex{2}()))
            marked = falses(size(m.nb)); marked[1,3] = true
            m1, parent = refine_with_parents(m, marked)
            @test nel(m1) > nel(m)
            @test parent[1:nel(m)] == 1:nel(m)
            @test nnodes(m1) == conforming_nnodes(m1)
            @test mesh_area(m1) ≈ mesh_area(m) rtol=1e-12
        end
        @test_throws DimensionMismatch refine(mshsquare(2), falses(3, 4))
        mp = mshsquare(2); set_bnd_periodic!(mp, (1,2), 1)
        @test_throws ErrorException uniref(mp)

        # Boundary layers: each pass halves the bottom row.
        m = mshsquare(4, p=2)
        m1, ix = bndlayer_refine_with_elements(m, 3, 2)
        @test nel(m1) == 24
        @test length(ix) == 12
        @test all(maximum(m1.x[m1.el[:,i],2]) <= 0.25 + 1e-12 for i in ix)
        @test minimum(maximum(m1.x[m1.el[:,i],2]) for i in 1:nel(m1)) ≈ 1/16
        @test nel(bndlayer_refine(m, (1, 3))) > nel(bndlayer_refine(m, 3))
        @test bndlayer_refine(m, 3, 2).x == m1.x
        # Thin layers far from the origin keep all their nodes.
        m = mshsquare(4, p=2); m.x[:,1] .+= 1000
        m1 = bndlayer_refine(m, 3, 16)
        @test nnodes(m1) == conforming_nnodes(m1)
        @test minimum(maximum(m1.x[m1.el[:,i],2]) for i in 1:nel(m1)) ≈ 0.25 / 2^16
        @test mesh_area(m1) ≈ 1 rtol=1e-12
        x, el = unique_mesh_nodes([1000 0; 1000 1e-9; 1000 1e-9; 1000 2e-9], [1 3; 2 4])
        @test size(x, 1) == 3 && el == [1 2; 2 3]

        # NACA mesh: the element counts agree with 3DG's mknaca1msh(1, 3).
        if Sys.which("gmsh") !== nothing
            naca = joinpath(@__DIR__, "data", "naca.geo")
            m  = set_lobatto_nodes(rungmsh2msh(naca; p=3, cmdadd="-v 0"))
            m1 = uniref(m)
            m2 = bndlayer_refine(m1, 1, 3)
            @test nel(m1) == 6144 && nel(m2) == 6564
            @test nnodes(m2) == conforming_nnodes(m2)
            @test mesh_area(m2) ≈ mesh_area(m) rtol=1e-10
            @test sort(unique(bndtag.(filter(isboundary, m2.nb)))) == [1, 2]
        end
    end

    @testset "Airfoil meshes" begin
        H = HighOrderMeshes
        # Coordinates: sample files, Lednicer format, NACA 4-digit formula.
        X = airfoil_coordinates(:naca0012)
        @test size(X) == (105, 2) && X[1,:] == X[end,:] == [1, 0] && X[53,:] == [0, 0]
        X = airfoil_coordinates(:rae2822)
        @test size(X) == (129, 2) && X[1,:] == X[end,:] == [1, 0] && X[65,:] == [0, 0]
        mktempdir() do dir
            fname = joinpath(dir, "lednicer.dat")
            open(fname, "w") do io
                println(io, "RAE 2822\n65. 65.\n")
                foreach(i -> println(io, X[i,1], " ", X[i,2]), [65:-1:1; 65:129])
            end
            @test airfoil_coordinates(fname) == X
        end
        @test_throws ErrorException airfoil_coordinates(:nosuchfoil)
        X = naca4("0012")
        @test size(X) == (401, 2) && X[1,:] == X[end,:]
        @test X[1,:] ≈ [1, 0] atol=1e-15
        @test X[1:201,2] ≈ -X[401:-1:201,2]
        @test maximum(X[:,2]) ≈ 0.06 rtol=1e-3               # 12% thickness
        X = naca4("2412")
        ix = argmin(abs.(X[1:201,1] .- 0.4))                 # max camber 2% at 40%
        @test (X[ix,2] + X[402-ix,2]) / 2 ≈ 0.02 atol=1e-3
        @test_throws ErrorException naca4("241")

        # Not-a-knot splines reproduce cubics; the stretching hits its end spacings.
        t  = (0:19) ./ 19 .+ 0.02 .* sin.(1:20)            # nonuniform knots
        sp = H._CubicSpline(t, [t.^3 .- t  2t.^2 .+ 1])
        for s in range(t[1], t[end], length=33)
            @test H._spline(sp, s) ≈ [s^3 - s, 2s^2 + 1] atol=1e-12
            @test H._dspline(sp, s) ≈ [3s^2 - 1, 4s] atol=1e-10
        end
        u = H._stretching(20, 0.01, 0.002)
        @test u[1] == 0 && u[end] ≈ 1 && all(diff(u) .> 0)
        @test u[2] - u[1] ≈ 0.01 && u[end] - u[end-1] ≈ 0.002
        # An open trailing edge is closed at the midpoint of the gap.
        X = naca4("0012"); X[1,2] += 0.002; X[end,2] -= 0.002
        g = H._airfoil_spline(X)
        @test H._spline(g.sp, 0.0) ≈ [1, 0] atol=1e-15
        @test H._spline(g.sp, g.T[end]) ≈ [1, 0] atol=1e-15

        if Sys.which("gmsh") !== nothing
            function mesh_area(m)
                ξ, w = quadrature(elgeom(m), 3*porder(m) + 2)
                J = interpolate(dshapefcns(m.fe, ξ), dg_nodes(m))
                dets = J[:,:,1,1] .* J[:,:,2,2] .- J[:,:,1,2] .* J[:,:,2,1]
                @assert minimum(dets) > 0
                sum(w .* dets)
            end
            function foil_area(X)
                g = H._airfoil_spline(X)
                P = reduce(vcat, (H._spline(g.sp, τ)' for τ in range(0, g.T[end], length=100001)))
                sum(P[i,1]*P[i+1,2] - P[i+1,1]*P[i,2] for i in 1:size(P,1)-1) / 2
            end
            fmap = facemap(Block{2}())
            for (foil, X, aoa) in ((:naca0012, airfoil_coordinates(:naca0012), 0),
                                   (:rae2822, airfoil_coordinates(:rae2822), 2),
                                   (naca4("2412"), naca4("2412"), 4))
                R, nlayers = 10, 3
                m = mshairfoil(foil; aoa, R, hmax=2, nfoil=24, nbndlayers=nlayers)
                @test elgeom(m) == Block{2}() && porder(m) == 3
                @test sort(unique(bndtag.(filter(isboundary, m.nb)))) == [1, 2]
                # The area is exact up to the wall geometry, and no element is inverted.
                @test mesh_area(m) ≈ 6R^2 - foil_area(X) atol=1e-6
                # Wall faces: 24 per surface, lengths hle and hte at the ends; the
                # layers are exact halvings of the band layers (tband/nband).
                cel = m.el[corner_nodes(m.fe), :]
                wall = [ (j, e) for e in axes(m.nb,2), j in 1:4 if isboundary(m.nb[j,e]) && bndtag(m.nb[j,e]) == 1 ]
                @test length(wall) == 48
                len(j, e) = hypot((m.x[cel[fmap[1,j],e],:] - m.x[cel[fmap[2,j],e],:])...)
                xmid(j, e) = sum(m.x[cel[fmap[:,j],e],1]) / 2
                @test minimum(len(f...) for f in wall if xmid(f...) > 0.9) ≈ 0.004 rtol=1e-4
                @test minimum(len(f...) for f in wall if xmid(f...) < 0.1) ≈ 0.004 rtol=1e-2
                heights = map(filter(f -> 0.45 < xmid(f...) < 0.55, wall)) do (j, e)
                    a, b = m.x[cel[fmap[1,j],e],:], m.x[cel[fmap[2,j],e],:]
                    n = [a[2] - b[2], b[1] - a[1]] / hypot((b - a)...)
                    maximum(abs((m.x[k,:] - a)' * n) for k in setdiff(cel[:,e], cel[fmap[:,j],e]))
                end
                @test heights ≈ fill(0.02 / 2 / 2^nlayers, length(heights)) rtol=1e-3
                # The layers pass through the TE: its 4 elements are 2 band and 2 wedge ones.
                iTE = findfirst(i -> m.x[i,:] ≈ X[1,:], axes(m.x, 1))
                @test count(==(iTE), cel) == 4
            end
            @test_throws ErrorException airfoil_geo(; aoa=40)
            @test_throws ErrorException airfoil_geo(; nfoil=25)
        end
    end

    @testset "Visualization data" begin
        for eg in (Block{2}(), Simplex{2}())
            m = ex1mesh(nref=1, eg=eg)
            elem_lines, int_lines, bnd_lines, elem_mid = viz_mesh(m)
            @test size(elem_mid) == (nel(m), 2)
            @test size(elem_lines, 2) == 2 && any(isnan, elem_lines)
            @test size(bnd_lines, 1) > 0
            for u in (ex1solution(m), ex1solution(m, dg=false))
                allx, allu, allel = viz_solution(m, u)
                @test size(allx, 1) == length(allu)
                @test size(allel, 1) == 3 && maximum(allel) <= size(allx, 1)
            end
        end
        m1 = set_degree(mshline(3), 2)
        allx, allu, allel = viz_solution(m1, m1.x)
        @test size(allel, 1) == 2
    end

    @testset "Solver node sets" begin
        # Gauss-Radau rules: left endpoint included, exact to degree 2n-2
        for n in 2:8
            x, w = gauss_radau_quadrature(n)
            @test length(x) == n && x[1] == -1 && all(diff(x) .> 0) && x[end] < 1
            for k in 0:2n-2
                @test w' * x.^k ≈ (iseven(k) ? 2/(k+1) : 0.0) atol=1e-12
            end
            x1, w1 = gauss_radau01_quadrature(n)
            @test x1[1] == 0 && sum(w1) ≈ 1
            @test w1' * x1.^(2n-2) ≈ 1/(2n-1)
        end
        @test gauss_radau_quadrature(2)[1] ≈ [-1, 1/3]
        @test eltype(gauss_radau01_nodes(4, Float32)) == Float32
        @test_throws ArgumentError gauss_radau_nodes(1)

        # Non-conforming lines are rejected for mesh elements, accepted with check=false
        @test_throws ArgumentError FiniteElement(Block{2}(), gauss_radau01_nodes(4))      # not symmetric
        @test_throws ArgumentError FiniteElement(Block{2}(), gauss_legendre01_nodes(4))   # no vertices
        @test_throws ArgumentError set_ref_nodes(mshsquare(2), gauss_legendre01_nodes(3))
        @test FiniteElement(Block{2}(), gauss_legendre01_nodes(4); check=false) isa FiniteElement

        # The DG-SEM workflow on affine block meshes: metric terms and exact transfer
        for (m, line) in ((set_degree(mshsquare(3), 3), gauss_legendre01_nodes(4)),
                          (set_lobatto_nodes(set_degree(mshcube(2, 2, 2), 2)), gauss_radau01_nodes(3)))
            G, p, D = elgeom(m), porder(m), dim(m)
            fe_sol  = FiniteElement(G, line; check=false)
            @test porder(fe_sol) == p && nnodes(fe_sol) == nnodes(m.fe)

            ξs = ref_nodes(fe_sol)
            xs = interpolate(shapefcns(m.fe, ξs),  dg_nodes(m))     # nξ × nel × D
            J  = interpolate(dshapefcns(m.fe, ξs), dg_nodes(m))     # nξ × nel × D × D
            h  = 1 / round(Int, nel(m)^(1/D))                       # axis-aligned elements of size h
            @test size(xs) == (nnodes(fe_sol), nel(m), D) && size(J) == (nnodes(fe_sol), nel(m), D, D)
            for i in 1:D, j in 1:D
                @test all(isapprox.(J[:,:,i,j], i == j ? h : 0, atol=1e-12))
            end

            # a degree-p polynomial sampled at solver nodes transfers exactly to the mesh nodes
            fx(X) = X[:,:,1].^p .- 2 .* X[:,:,2].^(p-1) .* X[:,:,1] .+ 3
            u_sol  = fx(xs)
            u_mesh = interpolate(fe_sol, m.fe, u_sol)
            @test u_mesh ≈ fx(dg_nodes(m)) atol=1e-10
            @test interpolate(m.fe, fe_sol, u_mesh) ≈ u_sol atol=1e-10
            @test size(interpolate(fe_sol, m.fe, cat(u_sol, 2u_sol, dims=3))) == (nnodes(m.fe), nel(m), 2)

            if D == 2
                # plotting straight from the solver layout agrees with the transferred field
                ax1, au1, ael1 = viz_solution(m, u_sol; fe=fe_sol)
                ax2, au2, ael2 = viz_solution(m, u_mesh)
                @test ax1 ≈ ax2 && au1 ≈ au2 && ael1 == ael2
                @test_throws ErrorException viz_solution(m, u_sol[1:end-1, :]; fe=fe_sol)
                @test_throws ArgumentError viz_solution(m, u_sol; fe=FiniteElement(Simplex{2}(), p))
            end
        end

        # Simplices: any unisolvent node set works with check=false; exact in reference space
        m  = ex1mesh(eg=Simplex{2}(), nref=1)
        p  = porder(m)
        ξs = Float64.(equispaced_nodes(Simplex{2}(), p)) .* 0.8 .+ 0.05
        @test_throws ArgumentError FiniteElement(Simplex{2}(), ξs)
        fe_sol = FiniteElement(Simplex{2}(), ξs; check=false)
        g(ξ) = ξ[:,1].^p .- ξ[:,2].^(p-1) .* ξ[:,1] .+ 1          # degree-p polynomial in ξ
        u_sol  = repeat(g(ξs), 1, nel(m))
        u_mesh = interpolate(fe_sol, m.fe, u_sol)
        @test u_mesh ≈ repeat(g(ref_nodes(m.fe)), 1, nel(m)) atol=1e-10
        @test_throws ArgumentError interpolate(fe_sol, FiniteElement(Block{2}(), p), u_sol)
    end

    @testset "Curved boundaries" begin
        fd(p) = sqrt(sum(p.^2)) - 1                       # unit circle
        fe(p) = (p[1]/2)^2 + p[2]^2 - 1                   # ellipse, not a distance function
        fsq(p) = -minimum((p[1], 1 - p[1], p[2], 1 - p[2]))   # unit square

        # Projection: the closest point for a distance function, on the zero set otherwise
        X = [2.0 0.0; 0.5 0.5; -0.3 0.1; 0.0 -3.0]
        Y = project_points(X, fd)
        @test all(abs.(norm.(eachrow(Y)) .- 1) .< 1e-14)
        @test Y ≈ X ./ norm.(eachrow(X)) atol=1e-8     # finite difference gradients
        @test all(abs.(fe.(eachrow(project_points(X, fe)))) .< 1e-13)
        @test project_points(zeros(0, 2), fd) == zeros(0, 2)
        @test project_points([1.0 0.0 0.0], p -> norm(p) - 1) == [1.0 0.0 0.0]

        # Arclength: exact on straight polylines, endpoints kept
        @test interp_arclength([0.0 0.0; 1.0 0.0; 1.0 3.0], [0, 0.25, 0.5, 1]) == [0 0; 1 0; 1 1; 1 3]
        t  = [0, 0.1, 0.15, 0.6, 1.0]
        s  = [0, 0.3, 0.5, 0.75, 1]
        a, b = [1.0 2.0], [-1.0 5.0]
        @test interp_arclength(a .+ t .* (b - a), s) ≈ a .+ s .* (b - a) atol=1e-15
        @test interp_arclength([1.0 2.0], [0.0, 0.5]) == [1.0 2.0; 1.0 2.0]

        # Blending matrices: vertex values by N, face values by B[j] on face j and zero on the others
        for fe1 in (FiniteElement(Simplex{2}(), 4), FiniteElement(Block{2}(), gauss_lobatto01_nodes(5)))
            N, Nf, B = HOM._boundary_blending(fe1)
            f2n = mkface2nodes(fe1)
            cn  = corner_nodes(fe1)
            fmap = facemap(elgeom(fe1))
            @test N[cn,:] ≈ I
            @test Nf ≈ shapefcns(Block{1}(), ref_nodes(fe1, 1))
            for j in axes(f2n,2)
                s1 = vec(ref_nodes(fe1, 1))
                φ  = randn(length(s1)) .* (0 .< s1 .< 1)   # face displacement without the vertex part
                for k in axes(f2n,2)
                    @test B[j][f2n[:,k],:] * φ ≈ (k == j ? φ : zero(φ)) atol=1e-13
                end
            end
        end

        # Straight boundaries on the zero level set do not move
        for m in (mshsquare(3; eg=Simplex{2}(), p=3), mshsquare(3; p=3), set_lobatto_nodes(mshsquare(2; p=5)))
            @test curve_boundary(m, fsq).x ≈ m.x atol=1e-14
        end

        # Hexagon refined to a disk: nodes on the circle at equal angles, conforming,
        # positive Jacobians, area converging, and only the elements at the boundary move
        x6 = [0 0; [cospi(k/3) for k in 0:5] [sinpi(k/3) for k in 0:5]]
        hexagon = HighOrderMesh(x6, [1 1 1 1 1 1; 2 3 4 5 6 7; 3 4 5 6 7 2])
        for p in 1:4
            errs = map((1, 2)) do nref
                m0 = set_degree(uniref(hexagon, nref), p)
                m  = curve_boundary(m0, fd)
                @test m.el == m0.el && m.nb == m0.nb
                for X in boundary_face_nodes(m, 1)
                    @test all(abs.(norm.(eachrow(X)) .- 1) .< 1e-13)
                    θ = atan.(X[:,2], X[:,1])
                    θ = mod.(θ .- θ[1] .+ π, 2π) .- π
                    @test θ ≈ θ[end] .* vec(ref_nodes(m.fe, 1)) atol=1e-10
                end
                bnd = boundary_nodes(m0)
                touching = unique(m0.el[:, [ any(in(bnd), m0.el[:,e]) for e in 1:nel(m0) ]])
                @test m.x[setdiff(1:nnodes(m), touching),:] == m0.x[setdiff(1:nnodes(m), touching),:]
                area, q = volume_and_quality(m)
                @test q > 0.6
                abs(area - π)
            end
            p > 1 && @test errs[1] / errs[2] > 2^(p+1) * 0.9
        end

        # Quads, curving only the arc of a quarter disk with straight edges
        m0 = set_degree(mshcircle(2; p=1, shape=:quarter), 3)
        m  = curve_boundary(m0, fd, 3)
        @test all(abs.(norm.(eachrow(reduce(vcat, boundary_face_nodes(m, 3)))) .- 1) .< 1e-13)
        @test boundary_face_nodes(m, 1) ≈ boundary_face_nodes(m0, 1) atol=1e-15
        @test boundary_face_nodes(m, 2) ≈ boundary_face_nodes(m0, 2) atol=1e-15
        @test volume_and_quality(m)[1] ≈ volume_and_quality(mshcircle(2; p=3, shape=:quarter))[1] rtol=1e-12

        # DistMesh: straight triangles in, curved high-order mesh out
        dm = DistMesh.distmesh2d(fd, DistMesh.huniform, 0.3, ((-1,-1), (1,1)))
        m1 = HighOrderMesh(dm)
        @test m1 isa HighOrderMesh{2,Simplex{2},Float64}
        @test nel(m1) == length(dm.t) && nnodes(m1) == length(dm.p)
        @test volume_and_quality(m1)[2] > 0
        m4 = curve_boundary(set_degree(m1, 4), fd)
        @test all(abs.(norm.(eachrow(m4.x[boundary_nodes(m4),:])) .- 1) .< 1e-13)
        area, q = volume_and_quality(m4)
        @test area ≈ π rtol=1e-6
        @test q > 0.5
    end
end
end
