using StatsBase, LinearAlgebra  

"""Generates n normal deviates with mean mu and standard
deviation sigma."""
function normal_deviate(mu::T, sigma::T, n::Int) where T <: Real
  x = zeros(T, n + 1)
  y = zeros(T, 2)
  for i = 1:2:n
    reject = true
    while reject
      y  = 2 * rand(T, 2) .- one(T)
      s = norm(y)^2
      if s < one(T) && s > zero(T)
        y = sqrt(-2 * log(s) / s) .* y
        reject = false
      end
    end
    x[i:i + 1] = y
  end
  x = sigma .* x .+ mu
  return x[1:n]
end

(mu, sigma, n) = (1.0, 2.0, 100000);
x = normal_deviate(mu, sigma, n);
describe(x)
