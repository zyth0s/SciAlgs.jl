"""Minimizes cos(x) by an MM algorithm and Newwton's method."""
function cosine_min(x)
  y = copy(x)
  for n = 0:6
    println(n,"  ",x,"  ",y)
    x = x + sin(x)  # MM update
    y = y - sin(y) / cos(y) # Newton update
  end
  return (x, y)
end

x = 2.0
(x, y) = cosine_min(x)
