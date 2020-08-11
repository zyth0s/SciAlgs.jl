
using SciAlgs.TheBook: babylonian

@testset "THE BOOK: Babylonian method" begin

   for tolexp in 1:7

      tol = exp10(-tolexp)

      @test  isapprox(babylonian(π^2,tol), π, atol=tol)
      @test  isapprox(babylonian(2.0^2,tol), 2, atol=tol)
   end
end

