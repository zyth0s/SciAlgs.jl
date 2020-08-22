# Author: Kenneth Lange @ University of California, Los Angeles 

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
    i, d = peek(pq) # retrieve the minimal remaining distance
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
  distance, predecessor
end
