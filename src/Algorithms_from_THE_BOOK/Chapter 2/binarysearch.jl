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
  return 0
end

x = ['a', 'b', 'd', 'f', 'g'];
println(binary_search(x, 'f'))
x = [1, 2, 4, 7, 9];
println(binary_search(x, 3))

