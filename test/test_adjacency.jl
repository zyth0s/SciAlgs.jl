
using SciAlgs.TheBook: adjacency_to_neighborhood

@testset "THE BOOK: Adjacency" begin

   #      1 2 3 4 5 6 7
   A = [[ 0 1 1 0 0 0 0]; # 1
        [ 1 0 1 0 0 0 0]; # 2
        [ 1 1 0 0 0 0 0]; # 3
        [ 0 0 0 0 1 1 0]; # 4
        [ 0 0 0 1 0 1 0]; # 5
        [ 0 0 0 1 1 0 0]; # 6
        [ 0 0 0 0 0 0 0]] # 7
   neighbor, weight = adjacency_to_neighborhood(A)
   @test [1,2] in neighbor
   @test [1,3] in neighbor
   @test [2,3] in neighbor
   @test [4,5] in neighbor
   @test [4,6] in neighbor
   @test [5,6] in neighbor
   #neighbor = [[2, 3], [1, 3], [1, 2], [5, 6], [4, 6], [4, 5], Int64[]]
   #weight   = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1], Int64[]]
end

