using Test

# -------------------------------------------------------------------
# 1. Core Library Tests
# -------------------------------------------------------------------
module TestCore
    using Test
    using HighOrderMeshes

    @testset verbose = true "Core HighOrderMeshes.jl" begin

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

    @testset "VTK Export" begin
        m = ex1mesh()
        u = ex1solution(m)
        
        mktempdir() do tmpdir
            filename = joinpath(tmpdir, "ex1.vtk")
            vtkwrite(filename, m, u)
            
            @test isfile(filename)
            @test filesize(filename) > 0
            header = open(readline, filename)
            @test startswith(header, "# vtk")
        end
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

        p = plot(m)
        check_plots(p)
        p = plot(m, u)
        check_plots(p)
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

        f = plot(m)
        check_makie(f)
        f = plot(m, u)
        check_makie(f)
    end
    end
else
    @info "Skipping Makie.jl tests. Set ENV[\"RUN_MAKIE_TESTS\"] = \"true\" to run them."
end
