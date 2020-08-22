
using SciAlgs.TheBook: adjacency_to_neighborhood, dijkstra

@testset "THE BOOK: Dijkstra's algorithm" begin

   #       1  2  3  4  5  6
   A = [[  0  7  9  0  0 14]; # 1
        [  7  0 10 15  0  0]; # 2
        [  9 10  0 11  0  2]; # 3
        [  0 15 11  0  6  0]; # 4
        [  0  0  0  6  0  9]; # 5
        [ 14  0  2  0  9  0]] # 6
   neighbor, weight = adjacency_to_neighborhood(A)

   distance, predecessor = dijkstra(neighbor, weight, 1) 
   @test distance ≈ [0.0, 7.0, 9.0, 20.0, 20.0, 11.0]
   @test predecessor == [0, 1, 1, 3, 6, 3]
end

