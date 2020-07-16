"""Collects neighborhoods and weights from an adjacency matrix A."""
function adjacency_to_neighborhood(A::AbstractMatrix)
#
  (nodes, T) = (size(A, 1), eltype(A))
  neighbor = [Vector{Int}() for i = 1:nodes]
  weight = [Vector{T}() for i = 1:nodes]
  for i = 1:nodes
    for j = 1:nodes
      if A[i, j] != zero(T)
        push!(neighbor[i], j)
        push!(weight[i], A[i, j])
      end
    end
  end
  return (neighbor, weight)
end

A = [[ 0 1 1 0 0 0 0]; [ 1 0 1 0 0 0 0]; [ 1 1 0 0 0 0 0];
[ 0 0 0 0 1 1 0]; [ 0 0 0 1 0 1 0]; [ 0 0 0 1 1 0 0];
[ 0 0 0 0 0 0 0]]; 
(neighbor, weight) = adjacency_to_neighborhood(A);

