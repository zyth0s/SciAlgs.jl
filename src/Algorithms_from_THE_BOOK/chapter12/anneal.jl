# Author: Kenneth Lange @ University of California, Los Angeles 

using Distances, Random

"""Solves traveling salesman problem by simulated annealing."""
function anneal(dist::Matrix{T}) where T <: Real
  n = size(dist, 1)
  steps = 50 * n^2 # annealing steps
  high_temp, low_temp = exp(8.0), exp(-6.5)
  rate = (low_temp / high_temp)^(1 / (steps - 1))
  temp = high_temp / rate
  path = randperm(n) # random initial path
  cost = dist[path[n], path[1]] # initialize cost
  for i = 2:n
    cost = cost + dist[path[i - 1], path[i]]
  end
  for i = 1:steps # commence annealing
    temp = temp * rate
    a, b = rand(2:n, 2)
    if a > b
      a, b = b, a
    elseif a == b
      continue
    end
    m = mod(b, length(path)) + 1 # compute change in cost
    c = dist[path[a - 1], path[b]] + dist[path[a], path[m]]
    c = c - dist[path[a - 1], path[a]] - dist[path[b], path[m]]
    if c < zero(T) || rand() < exp(- c / temp) # Metropolis
      cost = cost + c
      for i = 0:div(b - a, 2) # change the path
        path[a + i], path[b - i] = path[b - i], path[a + i]
      end
    end
  end
  cost, path
end

n = 50
X = rand(2, n)
dist = pairwise(Euclidean(), X, dims = 2)
cost, path = anneal(dist)
