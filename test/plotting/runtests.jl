# Run from the repository root with:
#
#   julia --project=test/plotting -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("test/plotting/runtests.jl")'

using Test
using HighOrderMeshes

@testset verbose = true "HighOrderMeshes plotting extensions" begin
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

        # plot! on an existing axis adds to it rather than replacing it
        f = Makie.Figure()
        ax = Makie.Axis(f[1,1])
        nplots0 = length(ax.scene.plots)
        Makie.plot!(ax, m)
        @test length(ax.scene.plots) > nplots0
        nplots1 = length(ax.scene.plots)
        Makie.plot!(ax, m, u)
        @test length(ax.scene.plots) > nplots1

        # plot(m, u) produces a Colorbar in 2D whose limits match the plotted
        # (interpolated) values, not just the nodal values
        f2 = Makie.plot(m, u)
        cbs = filter(x -> x isa Makie.Colorbar, f2.content)
        @test length(cbs) == 1
        allu = viz_solution(m, u)[2]
        @test collect(cbs[1].limits[]) ≈ collect(extrema(allu)) rtol=1e-6

        # a field on solver nodes plots through the fe attribute, and agrees with
        # the same field transferred to the mesh nodes
        fe_sol = FiniteElement(Block{2}(), gauss_legendre01_nodes(porder(m)+1); check=false)
        u_sol  = interpolate(m.fe, fe_sol, u)
        f3 = Makie.plot(m, u_sol; fe=fe_sol)
        check_makie(f3)
        f4 = Makie.plot(m, interpolate(fe_sol, m.fe, u_sol))
        lims(f) = collect(filter(x -> x isa Makie.Colorbar, f.content)[1].limits[])
        @test lims(f3) ≈ lims(f4) rtol=1e-6
    end
end
