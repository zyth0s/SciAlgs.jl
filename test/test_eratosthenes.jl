
using SciAlgs.TheBook: eratosthenes

@testset "THE BOOK: Eratosthenes Sieve" begin

      prime_list = eratosthenes(100)
      #@test   2 in prime_list
      #@test   3 in prime_list
      #@test   5 in prime_list
      #@test   7 in prime_list
      #@test  11 in prime_list
      #@test  13 in prime_list
      @test prime_list == [ 2, 3, 5, 7, 11, 13, 17, 19,
                           23, 29, 31, 37, 41, 43, 47,
                           53, 59, 61, 67, 71, 73, 79,
                           83, 89, 97] 
end

