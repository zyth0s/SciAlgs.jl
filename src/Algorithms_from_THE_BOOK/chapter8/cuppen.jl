# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra

""" Computes the roots of the quadratic a*x^2+b*x+c when |a|>0."""
function quadratic(a::T, b::T, c::T) where T <: Real
  d = b^2 - 4a * c # discriminant
  if d > zero(T) 
    if b >= zero(T)
      r1 = (-b - sqrt(d)) / (2a)
    else
      r1 = (-b + sqrt(d)) / (2a)
    end
    r2 = c / (r1 * a)
    return r1, r2
  else
    return (-b + sqrt(d + 0im)) / (2a), (-b - sqrt(d + 0im)) / (2a)
  end
end

"""Finds the root of the surrogate for the secular function."""
function approx_secular(d::Vector{T}, w::Vector{T}, mu::T,
  lambda::T, j::Int) where T <: Real
  n = length(d)
  f1, f2, c = zeros(T, 2), zeros(T, 2), zeros(T, 3)
  for i = 1:j # values and derivatives of the secular sums
    f1[1] = f1[1] + w[i] / (d[i] - lambda)
    f1[2] = f1[2] + w[i] / (d[i] - lambda)^2
  end
  for i = (j + 1):n
    f2[1] = f2[1] + w[i] / (d[i] - lambda)
    f2[2] = f2[2] + w[i] / (d[i] - lambda)^2
  end
  c[3] = mu # coefficients of the surrogate
  if j == 0
    c[2] = f2[2] * (d[j + 1] - lambda)^2
    c[3] = c[3] + f2[1] - f2[2] * (d[j + 1] - lambda)
    return (c[2] + c[3] * d[j + 1]) / c[3], mu + f1[1] + f2[1]
  elseif j == n
    c[1] = f1[2] * (d[j] - lambda)^2
    c[3] = c[3] + f1[1] - f1[2] * (d[j] - lambda)
    return (c[1] + c[3] * d[j]) / c[3], mu + f1[1] + f2[1]
  else
    c[1] = f1[2] * (d[j] - lambda)^2
    c[3] = c[3] + f1[1] - f1[2] * (d[j] - lambda)
    c[2] = f2[2] * (d[j + 1] - lambda)^2
    c[3] = c[3] + f2[1] - f2[2] * (d[j + 1] - lambda)
    p = c[3] # coefficients of the quadratic
    q = - (c[1] + c[2] + c[3] * (d[j] + d[j + 1]))
    r = c[1] * d[j + 1] + c[2] * d[j] + c[3] * d[j] *  d[j + 1]
    root1, root2 = quadratic(p, q, r) # solve the quadratic
    if Real(root1) > d[j] && Real(root1) < d[j + 1]
      return Real(root1), mu + f1[1] + f2[1]
    else
      return Real(root2), mu + f1[1] + f2[1]
    end
  end
end

"""Solves the secular equation by the LAPACK method."""
function solve_secular(d::Vector{T}, w::Vector{T}, rho::T, 
  tol::T) where T <: Real
  n = length(d)
  mu = one(T) / rho
  lambda, value = zeros(T, n), zeros(T, n)
  for i = 1:n
    if mu < zero(T) # choose the interval for ith root 
      j = i - 1
    else
      j = i
    end
    if j == 0 # initialize ith root and its bounds
      lambda[i], a, b = d[1] - one(T), -Inf, d[1] # on the left
    elseif j == n
      lambda[i], a, b = d[n] + one(T), d[n], Inf # on the right
    else
      lambda[i], a, b = (d[j] + d[j + 1]) / 2, d[j], d[j + 1] # midpoint
    end
    for k = 1:10 # iterate until convergence
      lambda[i], value[i] = approx_secular(d, w, mu, lambda[i], j)
      if minimum([abs(value[i]), lambda[i] - a, b - lambda[i]]) < tol
        break
      end
    end
  end
  lambda, value # roots and secular function values
end

"""Finds the eigenvector w associated with the eigenvalue alpha
of the matrix Diagonal(d) + r u * u'."""
function eigenvector(d::Vector{T}, u::Vector{T}, alpha::T, tol::T,
  i::Int) where T <: Real
  n, c = length(d), norm([d, u])
  w = zeros(T, n)
  if abs(u[i]) < tol * c # deflation cases first
    w[i] = one(T)
  elseif (i + 1 <= n) && (abs(d[i] - d[i + 1]) < tol * c)
    w[i] = one(T)
    w[i + 1] = - u[i] / u[i + 1]
  else # general case
    for j = 1:n
      w[j] = u[j] / (d[j] - alpha)
    end
  end
  w = w / norm(w)
end

"""Implements Cuppen's divide and conquer method of spectral 
decomposition for the symmetric tridiagonal matrix with 
diagonal a and subdiagonal b."""
function cuppen(a::Vector{T}, b::Vector{T}, tol::T) where T <: Real
  n = length(a)
  if n == 1 # trivial case
    return a[1] * ones(T, 1), ones(T, 1, 1)
  else
    V = zeros(T, n, n)
    k = div(n, 2) # divide before conquer
    rho = b[k] # multiplier of rank-1 perturbation
    (a1, b1) = ([a[1:k - 1]; a[k] - rho], b[1:k - 1])
    (a2, b2) = ([a[k + 1] - rho; a[k + 2:end]], b[k + 1:end])
    (e1, Q1) = cuppen(a1, b1, tol) # upper block decomposition
    (e2, Q2) = cuppen(a2, b2, tol) # lower block decomposition
    d = [e1; e2] # concatenate eigenvalues from blocks
    u = [vec(Q1[end, :]); vec(Q2[1, :])] # perturbation vector
    p = sortperm(d) # prepare to sort d 
    lambda, value = solve_secular(d[p], u[p].^2 , rho, tol)
    U = zeros(T, n, n)
    for i = 1:n # compute eigenvectors of the rank-1 perturbation
      U[:, i] = eigenvector(d, u, lambda[i], tol, i)
    end
    V[1:k, :] = Q1 * U[1:k, :] # convert eigenvectors
    V[k + 1:end,:] = Q2 * U[k + 1:end, :]
  end
  lambda, V
end

n, tol = 25, 1e-11
a, b = rand(n), randn(n - 1)
asave, bsave = copy(a), copy(b)
T = SymTridiagonal(a, b)
lambda, V = cuppen(a, b, tol)
println("norm test = ",norm(T * V - V * Diagonal(lambda)))
println("lambda = ",lambda)
println(eigvals(T))
