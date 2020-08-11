# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra

"""Extracts the Cholesky decomposition L of the symmetric positive
definite matrix A. The lower triangle of A is overwritten."""
function cholesky(A::Matrix{T}) where T <: Real
  n = size(A, 1)
  for k = 1:n
    A[k, k] = sqrt(A[k, k])
    for j = (k + 1):n
      A[j, k] = A[j, k] / A[k, k]
    end
    for i = (k + 1):n
      A[i:n, i] = A[i:n, i] - A[i, k] * A[i:n, k] 
    end
  end
  A
end

"""Solves the equation A*x = b using the Cholesky decomposition 
L of A. The vector b is overwritten."""
function choleskysolve(L::Matrix{T}, b::Vector{T}) where T <: Real
  n = size(L, 1)
  x = zeros(T, n)
  for j = 1:(n - 1) # forward solve
    b[j] = b[j] / L[j, j]
    b[j + 1:n] = b[j + 1:n] - b[j] * L[j + 1:n,j]
  end
  b[n] = b[n] / L[n, n]
  for i = n:-1:1 # backward solve
    x[i] = b[i]
    x[i] = x[i] - dot(L[i + 1:n,i], x[i + 1:n])
    x[i] = x[i] / L[i, i]
  end
  x
end

n = 100
M = randn(n, n)
A = M * M'
b = randn(n)
Asave, bsave = copy(A), copy(b)
L = cholesky(A)
x = choleskysolve(L, b)
norm(bsave - Asave * x)
