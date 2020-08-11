# Author: Kenneth Lange @ University of California, Los Angeles 

"""Minimizes cos(x) by an MM algorithm and Newton's method."""
function cosine_min(x)
  y = copy(x)
  for n = 0:6
    println(n,"  ",x,"  ",y)
    x = x + sin(x)  # MM update
    y = y - sin(y) / cos(y) # Newton update
  end
  x, y
end

x = 2.0
x, y = cosine_min(x)
