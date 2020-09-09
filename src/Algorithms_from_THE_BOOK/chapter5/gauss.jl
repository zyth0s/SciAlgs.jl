# Author: Kenneth Lange @ University of California, Los Angeles 

"""Overwrites the LU decomposition of A in A. Partial pivoting
is performed to enhance numerical stability."""
function LUdecomposition(A::Matrix{T}) where T <: Real
  n = size(A, 1)
  perm = collect(1:n) # identity permutation
  for k = 1:(n - 1)
    kp = argmax(abs.(A[k:n, k])) + k - 1 # find pivot row
    perm[k], perm[kp] = perm[kp], perm[k] # revise permutation
    for j = 1:n # interchange rows
      A[k, j], A[kp, j] = A[kp, j], A[k, j]
    end
    for i = (k + 1):n # compute multipliers
      A[i, k] = A[i, k] / A[k, k]
    end
    for j = (k + 1):n # adjust rows
      for i = (k + 1):n
        A[i, j] = A[i, j] - A[i, k] * A[k, j]
      end
    end
  end
  A, perm
end

"""Solves the equation A*x = b using the LU decomposition of A."""
function LUsolve(LU::Matrix{T}, perm::Vector{Int}, 
  b::Vector{T}) where T <: Real
  n = size(LU, 1)
  x = zeros(T, n)
  z = zeros(T, n)
  for i = 1:n # solve L*z = Perm*b
    z[i] = b[perm[i]]
    for j = 1:(i - 1)
      z[i] = z[i] - LU[i, j] * z[j]
    end
  end
  for i = n:-1:1 # solve U*x = z
    x[i] = z[i]
    for j = i + 1:n
      x[i] = x[i] - LU[i, j] * x[j]
    end
    x[i] = x[i] / LU[i, i]
  end
  x
end

"""Computes the inverse C of A. A is destroyed in the process."""
function inverse(A::Matrix{T}) where T <: Real
  n = size(A, 1)
  b = zeros(T, n)
  C = similar(A)
  LU, perm = LUdecomposition(A) # fetch LU decomposition
  for i = 1:n
    b[i] = one(T)
    C[:, i] = LUsolve(LU, perm, b) # solve A*x = e_i
    b[i] = zero(T)
  end
  C
end
