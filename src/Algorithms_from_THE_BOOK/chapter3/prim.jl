# Author: Kenneth Lange @ University of California, Los Angeles 

using DataStructures

"""Implements Prim's algorithm for finding a minimal
spanning tree of a graph."""
function prim(neighbor::Array{Array{Int, 1}, 1}, 
  weight::Array{Array{T, 1}, 1}) where T <: Number
#
  nodes = length(neighbor)
  visited = falses(nodes)
  mst_nodes, i = 1, 1 # initialize MST with node 1 as tip
  key = Array{Tuple{Int, Int}, 1}(undef, 0) # define keys
  priority = Array{Float64, 1}(undef, 0) # define priorities
  pq = PriorityQueue(zip(key, priority)) # initialize queue
  mst = Array{Tuple{Int, Int}, 1}(undef, 0) # minimum spanning tree
  while mst_nodes < nodes
    if !visited[i] # checked if the node has been visited
      visited[i] = true # mark the node as visited
      for k = 1:length(neighbor[i]) # add new edges to the queue
        j = neighbor[i][k]
        pq[(i, j)] = weight[i][k]
      end
    end
    (k, l), val = peek(pq) # choose lightest edge
    dequeue!(pq) # pop edge off queue
    if !visited[k] # if not part of tree, push the edge onto tree 
      push!(mst, (k, l))
      mst_nodes = mst_nodes + 1
      i = k # k is the new tip
    elseif !visited[l]
      push!(mst, (k, l))
      mst_nodes = mst_nodes + 1
      i = l # l is the new tip
    end
  end
  mst
end

"""Collects neighborhoods and weights from an adjacency matrix A."""
function adjacency_to_neighborhood(A::AbstractMatrix)
  nodes, T = size(A, 1), eltype(A)
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
  neighbor, weight
end

A = [0 2 0 6 0; 2 0 3 8 5; 0 3 0 0 7; 6 8 0 0 9; 0 5 7 9 0]
neighbor, weight = adjacency_to_neighborhood(A)
mst = prim(neighbor, weight)
