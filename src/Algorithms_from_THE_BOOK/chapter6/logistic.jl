# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra, Statistics

"""Performs logistic regression, where X is the design matrix, 
y is the response vector, n is the number of cases, and p is the 
number of parameters. """
function logistic(X::Matrix{T}, y::Vector{T}, tol::T) where T <:Real
  n, p = size(X)
  beta = zeros(T, p)
  beta[1] = log(mean(y) / (one(T) - mean(y))) # initialize intercept
  score, information = zeros(T, p), zeros(T, p, p)
  objective, old_objective = zero(T), zero(T)
  for iteration = 1:100
    fill!(score, zero(T))
    fill!(information, zero(T))
    objective = zero(T)
    for i = 1:n
      x = vec(X[i, :])
      inner = clamp(dot(x, beta), -20 * one(T), 20 * one(T))
      prob = exp(inner) / (exp(inner) + one(T))
      objective = objective + log(abs(prob + y[i] - one(T)))
      score = score + (y[i] - prob) * x
      information = information + prob * (one(T) - prob) * x * x'
    end
    old_objective = objective
    increment = information \ score  # scoring increment
    for step = 0:3 # step halving
      objective = zero(T)
      beta = beta + increment
      for i = 1:n
        x = vec(X[i, :])
        inner = clamp(dot(x, beta), -20 * one(T), 20 * one(T))
        prob = exp(inner) / (exp(inner) + one(T))
        objective = objective + log(abs(prob + y[i] - one(T)))
      end
      if old_objective < objective
        break
      else
        beta = beta - increment
        increment = increment / 2
      end
    end
    if iteration > 1 && abs(objective - old_objective) <= tol
      return (beta, iteration)
    end
  end
  beta, 100
end

n, p, tol = 100, 10, 1e-6
X, y, beta = randn(n, p), zeros(n), randn(p)
fill!(X[:, 1], 1.0) # intercept
for i = 1:n
  x = vec(X[i, :])
  inner = dot(x, beta)
  if rand() < exp(inner) / (exp(inner) + 1.0)
    y[i] = 1.0
  else
    y[i] = 0.0
  end
end
alpha = copy(beta)
beta, iteration = logistic(X, y, tol)
println(norm(beta - alpha))
