
using SciAlgs.TheBook: inverse

using Random: seed!

@testset "THE BOOK: LU decomposition" begin

   seed!(0)
   n = 100
   A = randn(n, n)
   Asave = copy(A)
   C = inverse(A)
   @test norm(Asave * C - I, 2) < 1e-12

end

