
using SciAlgs.TheBook: conjugategradients

using Random: seed!
using LinearAlgebra: norm, I
using SparseArrays: sprandn, opnorm

@testset "THE BOOK: Conjugate gradients" begin

   seed!(0)
   n, tol = 100, 1e-5
   A = sprandn(n, n, 0.01)
   rho = opnorm(A, 1)
   A = 0.5(A + A') + rho * I
   x = randn(n)
   b = A * x
   y = conjugategradients(A, b, tol)
   @test norm(x - y) < tol
end

