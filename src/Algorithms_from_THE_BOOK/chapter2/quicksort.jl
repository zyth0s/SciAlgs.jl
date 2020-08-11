# Author: Kenneth Lange @ University of California, Los Angeles 

"""Sorts the input list x from least to greatest by quicksort."""
function quicksort(x::Vector, left = 1, right = length(x))
  i = rand(left:right) # select a random splitting value
  split = x[i]
  x[left], x[i] = split, x[left]
  i = left
  for j = (left + 1):right # position the splitting value 
    if x[j] <= split
      i = i + 1
      x[i], x[j] = x[j], x[i]
    end
  end
  x[left], x[i] = x[i], split
  if i > left + 1 # sort to the left of the value
    quicksort(x, left, i - 1)
  end
  if i + 1 < right # sort to the right of the value
    quicksort(x, i + 1, right)
  end
end 

x = [5, 4, 3, 1, 2, 8, 7, 6, -1];
quicksort(x)
println(x)
x = ['a', 'c', 'd', 'b', 'f', 'e', 'h', 'g', 'y'];
quicksort(x)
println(x)
