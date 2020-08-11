# Author: Kenneth Lange @ University of California, Los Angeles 

using Distances

"""Implements kmeans clustering. The variable class should enter
with an initial guess of the classifications."""
function kmeans(X::Matrix{T}, class::Vector{Int}, 
  k::Integer) where T <: Real
#
  features, points = size(X)
  center, members = zeros(T, features, k), zeros(Int, k)
  switched = true
  while switched # iterate until memberships stabilize
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

k = 3
X = randn(100, 300)
X[:, 101:200] = X[:, 101:200] .+ 1.0
X[:, 201:300] = X[:, 201:300] .+ 2.0
class = rand(1:k, 300) # k classes randomly assigned
class, center = kmeans(X, class, k)
println(class)
