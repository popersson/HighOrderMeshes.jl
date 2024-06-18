using HighOrderMeshes
using Test

@testset "HighOrderMeshes.jl" begin
    for eg in (Block{2}(), Simplex{2}()), nref in 1:3
        @test porder(ex1mesh(nref=nref, eg=eg)) == 3
    end

    for eltpe in (Float64, Float32, Rational{Int})
        m = mshsquare(5, 3, eltpe)
        @test size(m.el) == (4,15)
        @test eltype(m.x) == eltpe
    end
end
