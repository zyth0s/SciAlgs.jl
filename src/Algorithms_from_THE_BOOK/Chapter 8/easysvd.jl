using LinearAlgebra

"""Extracts the svd of the matrix A."""
function easysvd(A::Matrix{T}) where T <: Real
  (m, n) = size(A)
  M = [zeros(n, n) A'; A zeros(m, m)]
  (D, V) = eigen(M)
  c = sqrt(2 * one(T))
  return (c * V[n + 1:end, 1:n] , D[1:n], c * V[1:n, 1:n])
end

(m, n) = (30, 20);
A = randn(m, n);
(U, S, V) = easysvd(A);
println(norm(A - U * Diagonal(S) * V'))

