# Author: Kenneth Lange @ University of California, Los Angeles 

""" Computes pi by Archimedes' algorithm."""
function archimedes(tol::T) where T <: Real
  a, b = 4 * one(T), 2 * sqrt(2 * one(T))
  while abs(a - b) > tol
    a = 2 * a * b / (a + b)
    b = sqrt(a * b)
  end
  a, b
end

upper, lower = archimedes(1e-6)
