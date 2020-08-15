
using SciAlgs.TheBook: heapsort!

@testset "THE BOOK: Heapsort algorithm" begin

   x = [5, 4, 3, 1, 2, 8, 7, 6, -1]
   x_sorted =  sort(x)
   heapsort!(x)
   @test x == x_sorted
   x = ['a', 'c', 'd', 'b', 'f', 'e', 'h', 'g', 'y']
   x_sorted =  sort(x)
   heapsort!(x)
   @test x == x_sorted
end

