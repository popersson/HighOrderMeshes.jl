using Test

# -------------------------------------------------------------------
# 1. Core Library Tests
# -------------------------------------------------------------------
module TestCore
using Test
using HighOrderMeshes
using LinearAlgebra
const HOM = HighOrderMeshes

# Evaluate the monomial ξ^α at every row of ξ.
monomial(ξ, α) = [ prod(ξ[k,d]^α[d] for d in eachindex(α); init=1.0) for k in axes(ξ,1) ]

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
                m       = gmshstr2msh(geo; porder=p, cmdadd="-v 0")
                fe      = m.fe
                xdg     = dg_nodes(m)
                corners = xdg[corner_nodes(fe), :, :]
                N1      = shapefcns(elgeom(m), ref_nodes(fe))
                pred    = interpolate(N1, corners)
                @test maximum(abs.(pred - xdg)) < 1e-10
            end
        end
    end

    @testset "CG Poisson (experimental)" begin
        # Solve -∇²u = 1 with zero Dirichlet boundary conditions on the unit circle
        for n = 1:4, porder = 1:4
            m = mshcircle(n, p=porder)
            pc = FEM_precomp(m)
            u,A,f = cg_poisson(m, pc, xy->1)
            uexact = (1 .- sum(m.x.^2,dims=2)) / 4
            error = maximum(abs.(u[:] - uexact[:]))
            # Assume O(h^{p+1}) convergence, with fitted constant (upper bound)
            error_bound = (0.2 / n) ^ (porder + 1)
            @test error < error_bound
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

        # NACA mesh: the element counts agree with 3DG's mknaca1msh(1, 3).
        if Sys.which("gmsh") !== nothing
            naca = joinpath(pkgdir(HighOrderMeshes), "examples", "naca", "naca.geo")
            m  = set_lobatto_nodes(rungmsh2msh(naca; porder=3, cmdadd="-v 0"))
            m1 = uniref(m)
            m2 = bndlayer_refine(m1, 1, 3)
            @test nel(m1) == 6144 && nel(m2) == 6564
            @test nnodes(m2) == conforming_nnodes(m2)
            @test mesh_area(m2) ≈ mesh_area(m) rtol=1e-10
            @test sort(unique(bndtag.(filter(isboundary, m2.nb)))) == [1, 2]
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
end
end
