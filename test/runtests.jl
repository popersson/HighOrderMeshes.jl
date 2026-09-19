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
        for (eg, pmax) in ((Simplex{2}(), 22), (Simplex{3}(), 15)), p in 1:pmax
            D = dim(eg)
            ξ, w = quadrature(eg, p)
            for α in multiindices(eg, p)
                exact = Float64(prod(factorial.(big.(α))) // factorial(big(sum(α) + D)))
                @test w' * monomial(ξ, α) ≈ exact atol=1e-11   # tabulated rules are accurate to ~1e-12
            end
        end
        @test eltype(quadrature(Block{2}(), 3, Float32)[2]) == Float32
        @test eltype(quadrature(Simplex{2}(), 3, Float32)[1]) == Float32

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
            m = gmsh2msh(fullname)
            @test elgeom(m) == eg
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

        @test_throws ErrorException uniref(set_degree(mshsquare(2), 2))
        @test nel(uniref(mshsquare(2), 2)) == 64
        @test length(boundary_nodes(set_degree(mshsquare(2), 3))) == 24
        @test length(boundary_nodes(set_degree(mshsquare(2), 3), 1)) == 7
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
end
end
