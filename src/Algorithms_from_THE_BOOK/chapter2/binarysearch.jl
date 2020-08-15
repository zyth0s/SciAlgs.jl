# Author: Kenneth Lange @ University of California, Los Angeles 

""" Conducts a binary search for the given value in the 
ordered list x."""
function binary_search(x::Vector, value)
  a = 1
  b = length(x)
  while a <= b
    m = div(a + b, 2)
    if x[m] > value 
      b = m - 1
    elseif x[m] < value 
      a = m + 1
    else 
      return m
    end
  end     
  0
end
