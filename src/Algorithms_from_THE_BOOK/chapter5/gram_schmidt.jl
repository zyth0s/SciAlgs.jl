# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra

"""Finds the QR decomposition of X. X should have more rows 
than columns."""
function gram_schmidt(X::Matrix{T}) where T <: Real
  n, p = size(X)
  R = zeros(T, p, p)
  Q = copy(X)
  for j = 1:p
    R[j, j] = norm(Q[:, j])
    Q[:, j] = Q[:, j] / R[j, j]
    for k = (j + 1):p
      R[j, k] = dot(Q[:, j], Q[:, k])
      Q[:, k] = Q[:, k] - Q[:, j] * R[j, k]
    end
  end
  Q, R
end

"""Solves the least squares problem of minimizing
||y - X * beta||^2 by the QR decomposition."""
function least_squares(X::Matrix{T}, y::Vector{T}) where T <: Real
  n, p = size(X)
  Q, R = gram_schmidt(X)
  beta = Q' * y
  for i = p:-1:1 # back substitution
    for j = (i + 1):p
      beta[i] = beta[i] - R[i, j] * beta[j]
    end
    beta[i] = beta[i] / R[i, i]
  end
  beta
end

n, p = 100, 20
X = randn(n, p)
y = randn(n)
beta = least_squares(X, y)
norm(beta - X \ y)
