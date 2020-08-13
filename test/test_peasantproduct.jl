
using SciAlgs.TheBook: peasantproduct

@testset "THE BOOK: Peasant product" begin

   @test  peasantproduct(10, 33) == 10 * 33
   a, b = rand.(Int, 2)
   @test  peasantproduct(a, b) == a * b
end

