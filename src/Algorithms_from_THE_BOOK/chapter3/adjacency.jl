# Author: Kenneth Lange @ University of California, Los Angeles 

"""Collects neighborhoods and weights from an adjacency matrix A."""
function adjacency_to_neighborhood(A::AbstractMatrix)

  nodes, T = size(A, 1), eltype(A)
  neighbor = [Vector{Int}() for i = 1:nodes]
  weight = [Vector{T}() for i = 1:nodes]
  for i = 1:nodes, j = 1:nodes
    if A[i, j] != zero(T)
      push!(neighbor[i], j)
      push!(weight[i], A[i, j])
    end
  end
  neighbor, weight
end
