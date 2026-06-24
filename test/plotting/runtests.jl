# Run from the repository root with:
#
#   julia --project=test/plotting -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("test/plotting/runtests.jl")'

using Test
using HighOrderMeshes

@testset verbose = true "HighOrderMeshes plotting extensions" begin
    @testset "Plots.jl Extension" begin
        import Plots
        import TriplotRecipes

        # Setup headless mode for GR (Plots backend)
        ENV["GKSwstype"] = "100"

        m = ex1mesh()
        u = ex1solution(m)

        function check_plots(p)
            @test p isa Plots.Plot
            mktempdir() do tmpdir
                path = joinpath(tmpdir, "test_plot.png")
                Plots.savefig(p, path)
                @test isfile(path)
                @test filesize(path) > 1000
            end
        end

        check_plots(Plots.plot(m, labels=:nodes))
        check_plots(Plots.plot(m, labels=:elements))
        check_plots(Plots.plot(m, u))
        check_plots(Plots.plot(m, u, mesh_edges=true))
        check_plots(Plots.plot(m, m.x, contours=10))

        m1 = set_degree(mshline(5), 3)
        check_plots(Plots.plot(m1, m1.x))
    end

    @testset "Makie Extension" begin
        import CairoMakie
        import Makie

        m = ex1mesh()
        u = ex1solution(m)

        function check_makie(f)
            @test f isa Makie.Figure
            ax = Makie.content(f[1,1]) # Get the axis from the figure layout
            @test !isempty(ax.scene.plots)
            mktempdir() do dir
                path = joinpath(dir, "test_makie.png")
                Makie.save(path, f)
                @test isfile(path)
            end
        end

        check_makie(Makie.plot(m, labels=:nodes))
        check_makie(Makie.plot(m, labels=:elements))
        check_makie(Makie.plot(m, u))
        check_makie(Makie.plot(m, u, mesh_edges=true))

        m1 = set_degree(mshline(5), 3)
        check_makie(Makie.plot(m1, m1.x))
    end
end
