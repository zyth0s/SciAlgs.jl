
using SciAlgs.TheBook: least_squares

using Random: seed!
using LinearAlgebra: norm

@testset "THE BOOK: Least squares by QR" begin

   seed!(0)
   n, p = 100, 20
   X = randn(n, p)
   y = randn(n)
   beta = least_squares(X, y)
   @test norm(beta - X \ y) < 3e-16
end

