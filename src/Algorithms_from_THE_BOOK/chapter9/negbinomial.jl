# Author: Kenneth Lange @ University of California, Los Angeles 

using Distributions

"""Estimates both parameters of the negative binomial 
distribution."""
function negbinomial(x::Vector{Int})
  m = length(x)
  p, r = 0.5, 1.0
  old_p, old_r = 0.5, 1.0
  avg = mean(x)
  for iteration = 1:500
    s = 0.0
    for i = 1:m
      for j = 0:(x[i] - 1)
        s = s + r / (r + j)
      end
    end
    r = -s / (m * log(p))
    p = r / (r + avg)
    if abs(p - old_p) + abs(r - old_r) < 1e-6
      break
    end
    (old_p, old_r) = (p, r)
  end
  p, r   
end

n, r, p = 100, 5.0, 0.25
x = rand(NegativeBinomial(r, p), n)
p, r = negbinomial(x)
obj = loglikelihood(NegativeBinomial(r, p), x)
println("r = ",r," p = ",p," obj = ",obj)
avg, ssq = (mean(x), var(x))
p, r = avg / ssq, avg^2 / (ssq - avg)
obj = loglikelihood(NegativeBinomial(r, p), x)
println("r = ",r," p = ",p," obj = ",obj)

