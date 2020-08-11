
using SciAlgs.TheBook: archimedes

@testset "THE BOOK: Archimedes method" begin

   for tolexp in 1:7

      tol = exp10(-tolexp)
      upper, lower = archimedes(tol)
      @test  lower < π < upper
   end
end

