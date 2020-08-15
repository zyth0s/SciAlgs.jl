
using SciAlgs.TheBook: quickselect!

@testset "THE BOOK: Quickselect method" begin

   k = 8
   x = [5, 4, 3, 1, 2, 8, 7, 6]
   xk = quickselect!(x, k)
   @test xk == 8

   k = 5
   x = ['a', 'c', 'd', 'b', 'f', 'e', 'h', 'g']
   xk = quickselect!(x, k)
   @test xk == 'e'
end

