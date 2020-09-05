
using SciAlgs: DNI, isvalid

@testset "DNI validation" begin

   dni = DNI(99999999, 'R')

   @test isvalid(dni)
end

