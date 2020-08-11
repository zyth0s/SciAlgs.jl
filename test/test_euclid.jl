
using SciAlgs.TheBook: euclid

@testset "THE BOOK: Euclidean algorithm" begin

   gcd_euclid = euclid(600, 220)
   @test gcd_euclid == gcd(600, 220)
end

