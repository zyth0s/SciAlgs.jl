using DataStructures

"""Implements Dykstra's algorithm for finding the shortest path
from a source node to every other node of a directed graph. 
The nodes should be numbered 1, 2,..."""
function dijkstra(neighbor::Array{Array{Int, 1}, 1}, 
  weight::Array{Array{T, 1}, 1}, source::Int) where T <: Number
#
  nodes = length(neighbor)
  node = collect(1:nodes)  # the nodes are numbered 1, 2, ...
  predecessor = zeros(Int, nodes)
  visited = falses(nodes)
  distance = zeros(nodes)
  fill!(distance, Inf)
  distance[source] = 0.0
  pq = PriorityQueue(zip(node, distance)) # priority queue
  while !isempty(pq)
    (i, d) = peek(pq) # retrieve the minimal remaining distance
    distance[i] = d
    visited[i] = true
    dequeue!(pq)  # pop the current node
    for k = 1:length(neighbor[i])
      j = neighbor[i][k]
      if !visited[j]
        dij = d + weight[i][k]
        if pq[j] > dij
          predecessor[j] = i
          pq[j] = dij  # adjust the provisional distance to j 
        end
      end
    end
  end
  return (distance, predecessor)
end

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

A = [[ 0 7 9 0 0 14]; [ 7 0 10 15 0 0]; [ 9 10 0 11 0 2];
[ 0 15 11 0 6 0]; [ 0 0 0 6 0 9]; [ 14 0 2 0 9 0]]; 
(neighbor, weight) = adjacency_to_neighborhood(A);
(distance, predecessor) = dijkstra(neighbor, weight, 1) 
 