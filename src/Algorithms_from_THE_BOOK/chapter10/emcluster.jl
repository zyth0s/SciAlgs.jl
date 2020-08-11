# Author: Kenneth Lange @ University of California, Los Angeles 

using Distances, Statistics, LinearAlgebra

"""Implements kmeans clustering. The variable class should enter
with an initial guess of the classifications."""
function kmeans(X::Matrix{T}, class::Vector{Int}, 
  k::Integer) where T <: Real
#
  features, points = size(X)
  center, members = zeros(T, features, k), zeros(Int, k)
  switched = true
  while switched # iterate until membership stabilizes
    fill!(center, zero(T)) # update centers
    fill!(members, 0) 
    for point = 1:points 
      i = class[point]
      center[:, i] = center[:, i] + X[:, point]
      members[i] = members[i] + 1
    end
    for j = 1:k
      center[:, j] = center[:, j] / max(members[j], 1)
    end
    switched = false # update classes
    dist = pairwise(Euclidean(), center, X, dims = 2) # fetch distances
    for point = 1:points
      j = argmin(dist[:, point]) # closest center
      if class[point] != j
        class[point] = j
        switched = true
      end
    end
  end 
  class, center
end
  
"""Performs EM cluster analysis. Assignment probabilities are
returned in the matrix assign. The data are stored in X, and
hard initial assignments enter in class."""
function emcluster(X::Matrix{T}, class::Vector{Int}, 
  tol::T) where T <: Real
#
  features, points = size(X) # initialize variables and arrays.
  classes = maximum(class)
  prior, center = zeros(T, classes), zeros(T, features, classes)
  assign = zeros(T, points, classes)
  for point = 1:points
    assign[point, class[point]] = one(T)
  end
  covar = cov(X')
  old_loglikelihood = -Inf
  for iteration = 1:1000 # enter the EM loop.
    center = X * assign # update centers and priors
    for i = 1:classes
      prior[i] = sum(assign[:, i])
      center[:, i] = center[:, i] / prior[i]
      prior[i] = prior[i] / points
    end
    fill!(covar, zero(T)) # update the covariance matrix
    for point = 1:points
      for i = 1:classes
        residual = X[:, point] - center[:, i]
        covar = covar + assign[point, i] * residual * residual'
      end
    end
    covar = covar / points
    inverse, lndet = inv(covar), logdet(covar)
    loglikelihood = zero(T) 
    for point = 1:points # update assignment probabilities
      for i = 1:classes
        residual = X[:, point] - center[:, i]
        quadratic = dot(residual, inverse * residual)
        assign[point, i] = prior[i] * exp(- (lndet + quadratic) / 2)
      end
      s = sum(assign[point, :]) # update the loglikelihood
      loglikelihood = loglikelihood + log(s)
      assign[point, :] = assign[point, :] / s
    end
    if loglikelihood < old_loglikelihood + tol # check convergence
      break
    end
    old_loglikelihood = loglikelihood
  end
  assign
end

k = 3
X = randn(100, 300)
X[:, 101:200] = X[:, 101:200] .+ 0.25
X[:, 201:300] = X[:, 201:300] .+ 0.50
guess = rand(1:k, 300)
class, center = kmeans(X, guess, k)
assign = emcluster(X, class, 1e-6)
for i = 1:300
  println(" kmeans class = ",class[i]," postprob = ", assign[i,:])
end
