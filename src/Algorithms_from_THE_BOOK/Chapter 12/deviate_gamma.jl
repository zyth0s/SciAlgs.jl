using Statistics

"""Generates n gamma deviates with shape parameter alpha
and intensity lambda."""
function gamma_deviate(alpha::T, lambda::T, n::Int) where T <: Real
#
  x = zeros(T, n)
  m = floor(Int, alpha)
  beta = alpha - m
  (y, z) = (zero(T), zero(T))
  for i = 1:n
    z = - log(prod(rand(T, m)))
    if beta <= one(T) / 10^6
      y = zero(T)
    elseif beta < one(T) / 10^3
      y = rand(T)^(one(T) / beta)
    else
      (r, s) = (one(T) / beta, beta - one(T))
      for trial = 1:1000
        u = rand(T, 2)
        y = - log(one(T) - u[1]^r)
        if u[2] <= (y / (one(T) - exp(- y)))^s
          exit
        end
      end
    end
    x[i] = (z + y) / lambda
  end
  return x
end

(alpha, lambda, n) = (3.5, 2.0, 10000);
x = gamma_deviate(alpha, lambda, n);
println(lambda * mean(x) / alpha, "  ", lambda^2 * std(x)^2 / alpha)
