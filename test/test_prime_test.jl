
using SciAlgs.TheBook: prime_test

using Primes

@testset "THE BOOK: Primality test" begin

   for i = 1:20
     n = rand(1:10000)
     prime = prime_test(n, 20)
     @test prime == isprime(n)
   end
end

