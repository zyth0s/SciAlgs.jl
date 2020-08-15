# Author: Kenneth Lange @ University of California, Los Angeles 

"""Finds the kth order statistic of a list x by the quick select 
method. The original x is destroyed in the process."""
function quickselect!(x::Vector, k::Int, left = 1, right = length(x))
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
  j = i - left + 1 # find the order statistic y
  if k == j
    y = x[i]
  elseif k < j
    y = quickselect!(x, k, left, i - 1)
  else
    y = quickselect!(x, k - j, i + 1, right)
  end
  y
end 
