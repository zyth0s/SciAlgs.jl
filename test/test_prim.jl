
using SciAlgs.TheBook: adjacency_to_neighborhood, prim

@testset "THE BOOK: Prim's algorithm" begin

   #    1 2 3 4 5
   A = [0 2 0 6 0; # 1
        2 0 3 8 5; # 2
        0 3 0 0 7; # 3
        6 8 0 0 9; # 4
        0 5 7 9 0] # 5
   neighbor, weight = adjacency_to_neighborhood(A)

   mst = prim(neighbor, weight)
   @test mst == [ (1, 2), (2, 3), (2, 5), (1, 4) ]

end

