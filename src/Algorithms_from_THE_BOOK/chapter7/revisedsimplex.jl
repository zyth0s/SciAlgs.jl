# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra, SparseArrays

"""Minimizes c'x subject to Ax = b and x >= 0. B enters
with the indices of a feasible vertex. The linear program
is assumed to be nondegenerate.""" 
function revised_simplex(A::AbstractMatrix{T}, b::Vector{T},  
  c::Vector{T}, B::Vector{Int}, tol::T) where T <: Real
  m, n = size(A)
  N = setdiff(collect(1:n), B)
  ABinv = inv(convert(Matrix{T}, A[:, B]))
  xB = ABinv * b
  for iteration = 1:10*m
    mu = c[N]' - (c[B]' * ABinv) * A[:, N]
    k = argmin(mu')
    if mu[k] > -tol # test for convergence
      return (dot(c[B], xB), B, xB)
    else
      d = ABinv * A[:, N[k]] # possible update directions
      p = findall(d .> tol)
      if isempty(p)
        return -Inf, nothing, nothing # unbounded below
      end
      t, i = findmin(xB[p] ./ d[p])
      B[p[i]], N[k] = N[k], B[p[i]] # revise vertex set
      xB = xB - t * d
      xB[p[i]] = t
      v = ABinv[p[i], :] / d[p[i]] # Sherman-Morrison update
      d[p[i]] = d[p[i]] - one(T)
      ABinv = ABinv - d * v'
    end
  end
end

"""Orchestrates the revised simplex method."""
function simplex_program(A::AbstractMatrix{T}, b::Vector{T},  
  c::Vector{T}, tol::T) where T <: Real
  m, n = size(A)
  for i = 1:m
    if b[i] < zero(T)
      A[i, :] = - A[i, :]
      b[i] = - b[i]
    end
  end
  A1 = [A I]
  c1 = [zeros(T, n); ones(T, m)]
  B = collect(n + 1:m + n)
  cost, B, xB = revised_simplex(A1, b, c1, B, tol) # first phase
  if cost > tol
    return "infeasible", Inf, cost, B, xB
  else
    cost, B, xB = revised_simplex(A, b, c, B, tol) # second phase
    if cost < -1e6
      return "unbounded", -Inf, B, xB
    else 
      return "solvable", cost, B, xB
    end
  end
end

# c = [-3.0; -2; 0; 0; 0]; # Wikipedia example
# A = [ 1.0 1 1 0 0; 2 1 0 1 0; 1 0 0 0 1 ]
# b = [ 40.0; 60; 30]
c = [-2.0; -3; -4; 0; 0] # Wikipedia example
A = [ 3.0 2 1 1 0; 2 5 3 0 1 ]
b = [ 10.0; 15 ]
A = sparse(A)
tol = 1e-5
status, cost, B, xB = simplex_program(A, b, c, tol)
