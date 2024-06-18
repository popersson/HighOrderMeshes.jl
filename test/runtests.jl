using HighOrderMeshes
using Test

include("../examples/sample-meshes.jl")

@testset "HighOrderMeshes.jl" begin
    @test porder(test_quad3()) == 3
end
