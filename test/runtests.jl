using CompariMotif
using Test
using Aqua
using JET

@testset "CompariMotif.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CompariMotif)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(CompariMotif; target_defined_modules = true)
    end
    # Write your tests here.
end
