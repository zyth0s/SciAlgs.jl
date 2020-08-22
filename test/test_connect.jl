
using SciAlgs.TheBook: connect

@testset "THE BOOK: Connect" begin

   #      1 2 3 4 5 6 7
   #A = [[ 0 1 1 0 0 0 0]; # 1
   #     [ 1 0 1 0 0 0 0]; # 2
   #     [ 1 1 0 0 0 0 0]; # 3
   #     [ 0 0 0 0 1 1 0]; # 4
   #     [ 0 0 0 1 0 1 0]; # 5
   #     [ 0 0 0 1 1 0 0]; # 6
   #     [ 0 0 0 0 0 0 0]] # 7
   #neighbor, weight = adjacency_to_neighborhood(A)
   neighbor = [[2, 3], [1, 3], [1, 2], [5, 6], [4, 6], [4, 5], Int64[]]
   weight   = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1], Int64[]]
   component, components = connect(neighbor)
   @test component == [1,  1,  1,  2,  2,  2,  3]
   @test components == 3
end

