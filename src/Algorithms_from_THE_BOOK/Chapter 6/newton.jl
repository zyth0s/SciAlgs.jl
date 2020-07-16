""" Implements Newton's method in one dimension for finding a 
root of the equation f(x) = 0. fp is the derivative of f(x)."""
function newton(f::Function, fp::Function, x::Real, tol::Real)
  (xold, xnew) = (x, zero(x))
  for iteration = 1:100
    xnew = xold - f(xold) / fp(xold)
    if abs(f(xnew)) < tol
      return (xnew, iteration)
    end
    xold = xnew
  end
  return (xnew, 100)
end

f(x) = x^3 - 2x - 5.0
fp(x) = 3x^2 - 2.0
(x, iterations) = newton(f, fp, 2.5, 1e-14)
f(x)
