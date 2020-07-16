using LinearAlgebra, SparseArrays

"""Implements the pagerank algorithm. Q is the transpose of the transition
probability matrix, u is the initial value of the equilibrium distribution, 
and 1-alpha is the probability of moving at each epoch according to Q."""
function pagerank(Q::AbstractMatrix{T}, u::Vector{T}, alpha::T,
  n::Int) where T <: Real
#
  c = alpha / length(u)
  for i = 1:n
    u = (one(T) - alpha) * (Q * u) .+ c
  end
  return u
end

"""Prepares data for Pagerank."""
function transpose_transition_matrix(m::Int)
#
  Q = spzeros(m, m)
  for i = 1:m
    neighbor = rand(1:m, 5) # each node has about 5 neighbors
    neighbor = setdiff(unique(neighbor), [i])
    for j in neighbor
      Q[j, i] = 1.0 / length(neighbor)
    end 
  end
  return Q
end

m = 100;
Q = transpose_transition_matrix(m);
u = rand(m);
u = u / sum(u);
(alpha, n) = (0.05, 20);
v = pagerank(Q, u, alpha, n);
R = (1.0 - alpha) * Q + (alpha / m) * ones(m, m);
norm(R * v - v)
