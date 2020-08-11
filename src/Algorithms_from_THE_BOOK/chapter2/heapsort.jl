# Author: Kenneth Lange @ University of California, Los Angeles 

"""Sorts the input list x from least to greatest by heapsort."""
function heapsort(x::Vector)
  n = length(x)
  for parent = div(n, 2):-1:1 # form the heap
    siftdown(x, parent, n)
  end
  for bottom = n:-1:2 # sort the heap
    (x[1], x[bottom]) = (x[bottom], x[1])
    siftdown(x, 1, bottom - 1)
  end
end
 
"""Moves an item down in the binary tree until it no longer
conflicts with the heap requirement."""
function siftdown(x::Vector, parent::Int, bottom::Int)
  parent_value = x[parent] 
  child = 2 * parent 
  while child <= bottom  
    if child < bottom && x[child] < x[child + 1]
      child = child + 1
    end 
    if x[child] <= parent_value
      break
    else
      x[div(child, 2)] = x[child]
      child = 2 * child
    end
  end
  x[div(child, 2)] = parent_value
  nothing
end

x = [5, 4, 3, 1, 2, 8, 7, 6, -1]
heapsort(x)
println(x)
x = ['a', 'c', 'd', 'b', 'f', 'e', 'h', 'g', 'y']
heapsort(x)
println(x)

