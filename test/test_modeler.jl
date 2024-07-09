@testset "modeler.jl" begin
    # simple test of convergence
    edm = rand_edm(10, 2)
    @test solve(edm, 4, "random", 2, GLPK) > 0

    # test exactness

end
