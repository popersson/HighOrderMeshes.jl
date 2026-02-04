using Test

# -------------------------------------------------------------------
# 1. Core Library Tests
# -------------------------------------------------------------------
module TestCore
using Test
using HighOrderMeshes

@testset verbose = true "Core HighOrderMeshes.jl" begin
    
    @testset "Element Geometries" begin
        for D in 1:3, G in [Simplex, Block]
            geom = G{D}()
                
            @testset "$G{$D}" begin
                # 1. Basic properties
                @test dim(geom) == D
                @test nvertices(geom) > 0
                @test nfaces(geom) > 0
                @test nedges(geom) >= 0
                
                # 2. Consistency Checks
                if D < 4
                    fmap = facemap(geom)
                    @test size(fmap,2) == nfaces(geom)
                    emap = edgemap(geom)
                    @test size(emap,2) == nedges(geom)
                end
                
                # 3. Inverse lookup test
                if !(D == 1 && G == Simplex) # Cannot tell Simplex from Block in 1D
                    nv = nvertices(geom)
                    @test HighOrderMeshes.find_elgeom(D, nv) == geom
                end
                
                # 4. Helper functions
                @test subgeom(geom, D-1) isa ElementGeometry{D-1}
                @test length(HighOrderMeshes.midpoint(geom)) == D
                
                # 5. Show methods, check that it prints to a buffer
                buf = IOBuffer()
                show(buf, geom)
                s = String(take!(buf))
                @test startswith(s, "ElementGeometry:")
                @test contains(s, "$(D)D")
            end
        end
        
        # Corner Cases / Error Handling
        @test_throws ErrorException HighOrderMeshes.find_elgeom(2, 100)
        @test_throws ErrorException HighOrderMeshes.plot_face_order(Simplex{4}())
        @test_throws ErrorException facemap(Simplex{4}())
        @test_throws ErrorException edgemap(Simplex{4}())
    end

    @testset "Basic ex1mesh properties" begin
        for eg in (Block{2}(), Simplex{2}()), nref in 1:3
            @test porder(ex1mesh(nref=nref, eg=eg)) == 3
        end
    end
        
    @testset "Basic mshsquare properties" begin
        for xtype in (Float64, Float32, Rational{Int})
            m = mshsquare(5, 3; T=xtype)
            @test size(m.el) == (4,15)
            @test eltype(m.x) == xtype
        end
    end
        
    @testset "IO - .hom and .txt formats" begin
        for eg in (Block{2}(), Simplex{2}()), nref in 0:1
            m = ex1mesh(nref=nref, eg=eg)
    
            mktempdir() do tmpdir
                tmp_path = joinpath(tmpdir, "test_mesh.hom")
                
                savemesh(tmp_path, m)
                m2 = loadmesh(tmp_path)
                
                @test m.x == m2.x && m.el == m2.el && m.nb == m2.nb &&
                    m.fe.ref_nodes[1] == m2.fe.ref_nodes[1]
            end
            
            mktempdir() do tmpdir
                m = set_degree(m, 1) # txt-format only supports p=1 for now
                tmp_path = joinpath(tmpdir, "test_mesh.txt")
                
                savemeshtxt(tmp_path, m)
                m2 = loadmeshtxt(tmp_path)
                
                @test m.x == m2.x && m.el == m2.el
            end
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
        rootdir = pkgdir(HighOrderMeshes)
        include(joinpath(rootdir, "examples", "fem", "assemble_utils.jl"))
        # Solve and plot -∇²u = 1 with zero Dirichlet boundary conditions on the unit circle
        for n = 1:4, porder = 1:4
            m = mshcircle(n, p=porder)
            pc = FEM_precomp(m)
            u,A,f = cg_poisson(m, pc, xy->1)
            uexact = (1 .- sum(m.x.^2,dims=2)) / 4
            error = maximum(abs.(u[:] - uexact[:]))
            # @show (n,porder,error)
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
        end
    end

    @testset "Mesh utilities" begin
        msh = mshsquare(5)
        set_bnd_periodic!(msh, (1,2), 1)    # Periodic left/right (x-direction)
        set_bnd_periodic!(msh, (3,4), 2)    # Periodic bottom/top (y-direction)
        @test minimum(first.(msh.nb)[:]) == 1  # No actual boundaries

        msh = ex1mesh(nref=1, eg=Block{2}())
        align_with_ldgswitch!(msh)
        sw = mkldgswitch(msh)
        @test all(sw[1,:] .== 1 .&& sw[3,:] .== 1)
    end

    @testset "Polynomials" begin
        # Quadrature
        maxdegree = (30,22,15)
        for D = 1:3
            for porder = 1:maxdegree[D]
                for (eg,exact) in ((Simplex{D}(), 1/(D+1)),
                                   (Block{D}(), 0.5))
                    gx,gw = quadrature(eg, porder)
                    # Trivial test int(x)
                    @test gw' * gx[:,1] ≈ exact
                end
            end
        end

        
        # Compute integral(x^4, x=[-1,1]) = 2/5
        # n = 3 => DoP = 5 => Exact
        x,w = gauss_legendre_quadrature(3)
        @test w' * x.^4 ≈ 2/5

        x,w = gauss_lobatto_quadrature(4)
        @test w' * x.^4 ≈ 2/5
        
        # Compute integral(x^4, x=[0,1]) = 1/5
        # n = 3 => DoP = 5 => Exact
        x,w = gauss_legendre01_quadrature(3)
        @test w' * x.^4 ≈ 1/5

        x,w = gauss_lobatto01_quadrature(4)
        @test w' * x.^4 ≈ 1/5
        
    end
        
end
end

# -------------------------------------------------------------------
# 2. Test Plots.jl Extension
# -------------------------------------------------------------------
if get(ENV, "RUN_PLOTS_TESTS", "false") == "true"
    @eval module TestPlots
    
    using Test
    using HighOrderMeshes
    using Plots
    using TriplotRecipes 

    # Setup headless mode for GR (Plots backend)
    ENV["GKSwstype"] = "100"

    @testset "Plots.jl Extension" begin
        m = ex1mesh()
        u = ex1solution(m)

        function check_plots(p)
            @test p isa Plots.Plot
            mktempdir() do tmpdir
                path = joinpath(tmpdir, "test_plot.png")
                savefig(p, path)
                @test isfile(path)
                @test filesize(path) > 1000
            end
        end

        f = plot(m, labels=:nodes)
        check_plots(f)
        f = plot(m, labels=:elements)
        check_plots(f)
        f = plot(m, u)
        check_plots(f)
        
        f = plot(m, u, mesh_edges=true)
        check_plots(f)

        f = plot(m, m.x, contours=10)
        check_plots(f)

        m1 = mshline(5)
        m1 = set_degree(m1, 3)
        f = plot(m1, m1.x)
        f = check_plots(f)
    end
    end
else
    @info "Skipping Plots.jl tests. Set ENV[\"RUN_PLOTS_TESTS\"] = \"true\" to run them."
end

# -------------------------------------------------------------------
# 3. Test Makie Extension
# -------------------------------------------------------------------
if get(ENV, "RUN_MAKIE_TESTS", "false") == "true"
    @eval module TestMakie
    using Test
    using HighOrderMeshes
    using CairoMakie 

    @testset "Makie Extension" begin
        m = ex1mesh()
        u = ex1solution(m)
        
        function check_makie(f)
            @test f isa Makie.Figure
            ax = content(f[1,1]) # Get the axis from the figure layout
            @test !isempty(ax.scene.plots) 
            mktempdir() do dir
                save(joinpath(dir, "test_makie.png"), f)
                @test isfile(joinpath(dir, "test_makie.png"))
            end
        end

        f = plot(m, labels=:nodes)
        check_makie(f)
        f = plot(m, labels=:elements)
        check_makie(f)
        f = plot(m, u)
        check_makie(f)
        
        f = plot(m, u, mesh_edges=true)
        check_makie(f)

        m1 = mshline(5)
        m1 = set_degree(m1, 3)
        f = plot(m1, m1.x)
        f = check_makie(f)
    end
    end
else
    @info "Skipping Makie.jl tests. Set ENV[\"RUN_MAKIE_TESTS\"] = \"true\" to run them."
end
