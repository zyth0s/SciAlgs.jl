using SciAlgs
using Test

import SciAlgs: my_f

my_f(2,1)

@testset "SciAlgs.jl" begin
   @test my_f(2,1) == 5
end

include("test_gauss_legendre.jl")
