""" Computes the square root of c >= 0."""
function babylonian(c::T, tol::T) where T <: Real
  x = one(T)  # start x at 1
  while abs(x^2 - c) > tol  # convergence test 
    x = (x + c / x) / 2
  end
  return x
end

root = babylonian(pi^2, 1e-10)

