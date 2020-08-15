# Author: Kenneth Lange @ University of California, Los Angeles 

""" Implements the bisection algorithm for finding a root of
the equation f(x)=0. The constants a < b should bracket a 0."""
function bisect(f::Function, a::T, b::T, tol::T) where T <: Real
  fa, fb = f(a), f(b)
  @assert(a < b && fa * fb <= zero(T)) # check for input error
  for iteration = 1:100
    m = (a + b) / 2
    fm = f(m)
    if abs(fm) < tol
      return m, iteration
    end
    if fa * fm < zero(T) 
      b, fb = m, fm
    else
      a, fa = m, fm
    end
  end
  (a + b) / 2, 100
end

f(x) = x^3 - 5x + 1.0
x, iteration = bisect(f, 0.0, 2.0, 1e-14)
