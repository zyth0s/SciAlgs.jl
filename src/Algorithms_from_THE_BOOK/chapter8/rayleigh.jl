# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra

"""Finds an eigenvalue and corresponding eigenvector close
to the initial vector v."""
function rayleigh(A::Matrix{T}, v::Vector{T}, tol::T) where T <: Real
  v = v / norm(v)
  mu = dot(v, A * v)
  for i = 1:10
    u = (A - mu * I) \ v
    u = u / norm(u)
    mu = dot(u, A * u)
    test = norm(u - v)
    v = copy(u)
    if test < tol
      break
    end  
  end
  mu, v
end

"""Delivers the dominant eigenvector of the matrix A."""
function power_method(A::Matrix{T}) where T <: Real
  v = randn(n)
  for i = 1:100 
    v = A * v
    v = v / norm(v)
  end
  v
end

n, tol = 100, 1e-10
A = randn(n, n)
A = A + A'
v = power_method(A)
mu, v = rayleigh(A, v, tol) # improve the power estimate
println("error = ",norm(A * v - mu * v)," mu = ",mu)
