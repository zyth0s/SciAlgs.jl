
using SciAlgs.TheBook: cholesky, choleskysolve

using Random: seed!
using LinearAlgebra: norm

@testset "THE BOOK: Cholesky decomposition" begin

   seed!(0)
   n = 100
   M = randn(n, n)
   A = M * M'
   b = randn(n)
   Asave, bsave = copy(A), copy(b)
   L = cholesky(A)
   x = choleskysolve(L, b)
   @test norm(bsave - Asave * x) < 1e-11
end

