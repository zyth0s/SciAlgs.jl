# Author: Kenneth Lange @ University of California, Los Angeles 

using Roots, Distributions, LinearAlgebra, SpecialFunctions

"""Finds the maximum likelihood estimates of the parameters 
of the Dirichlet distribution."""
function dirichlet(X::Matrix{T}) where T <: Real
  p = size(X, 1)
  avglog = mean(log.(X), dims = 2)
  lambda, old_lambda = ones(p), ones(p)
  for iteration = 1:500
    c = digamma(sum(lambda))
    for i = 1:p
      f(r) = avglog[i] + c - digamma(r)
      lambda[i] = find_zero(f, 1.0) # solve the 1-d update eqs
    end
    if norm(lambda - old_lambda) < 1e-6
      break
    end
    old_lambda = copy(lambda)
  end
  return lambda
end 

n, p = 1000, 3
lambda = 5rand(p)
X = rand(Dirichlet(lambda), n)
lambda = fit_mle(Dirichlet, X)
println("Package estimated parameters = ",lambda)
lambda = dirichlet(X)
println("Direct estimated parameters = ",lambda)
