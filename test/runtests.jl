using HighOrderMeshes
using Test

@testset "HighOrderMeshes.jl" begin
    for eg in (Block{2}(), Simplex{2}()), nref in 1:3
        @test porder(ex1mesh(nref=nref, eg=eg)) == 3
    end

    for eltpe in (Float64, Float32, Rational{Int})
        m = mshsquare(5, 3; T=eltpe)
        @test size(m.el) == (4,15)
        @test eltype(m.x) == eltpe
    end

    rootdir = pkgdir(HighOrderMeshes)
    for (filename,eg) in (("circle_tris.msh", Simplex{2}()),
                          ("square_tris.msh", Simplex{2}()),
                          ("circle_quads.msh", Block{2}()),
                          ("square_quads.msh", Block{2}()))
        fullname = joinpath(rootdir, "examples/gmsh", filename)
        m = gmsh2msh(fullname)
        @test elgeom(m) == eg
    end
end
