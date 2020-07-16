"""Calculates calculates log(x) with a cubic rate of convergence."""
function fastlog(x::Real)
  y = exponent(x) * 0.6931471805599453
  for n = 1:3
    z = exp(y)
    y = y - 2 * (z - x) / (z + x)
  end
  return y
end

fastlog(1e20 * pi) - log(1e20 * pi)
