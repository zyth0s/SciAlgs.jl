
using SciAlgs.TheBook: binary_search

@testset "THE BOOK: Binary search" begin

   x = ['a', 'b', 'd', 'f', 'g']
   @test binary_search(x, 'f') == 4

   x = [1, 2, 4, 7, 9]
   @test binary_search(x, 3) == 0
end

