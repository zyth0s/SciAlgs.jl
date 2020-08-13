
using SciAlgs.TheBook: quadratic

@testset "THE BOOK: Roots of quadratic equation" begin
   # x² - 2x + 1 = 0; (x - 1)(x - 1) = 0
   a, b, c = 1.0, -2.0, 1.0
   r1, r2 = quadratic(a, b, c)
   @test  real(r1) ≈ 1
   @test  imag(r1) ≈ 0
   @test  real(r2) ≈ 1
   @test  imag(r2) ≈ 0
end

