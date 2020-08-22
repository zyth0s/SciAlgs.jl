# Author: Kenneth Lange @ University of California, Los Angeles 

"""Finds the connected components of a graph.  Nodes should be 
numbered 1,2,..."""
function connect(neighbor::Array{Array{Int, 1}, 1})

  nodes = length(neighbor)
  component = zeros(Int, nodes)
  components = 0
  for i = 1:nodes
    if component[i] > 0 continue end
    components = components + 1
    component[i] = components
    visit!(neighbor, component, i)
  end
  component, components
end

"""Recursively assigns components by depth first search."""
function visit!(neighbor::Array{Array{Int,1},1},
                component::Vector{Int}, i::Int)

  for j in neighbor[i]
    if component[j] > 0 continue end
    component[j] = component[i]
    visit!(neighbor, component, j)
  end
end
