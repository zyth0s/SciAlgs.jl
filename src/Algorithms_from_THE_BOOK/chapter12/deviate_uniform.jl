# Author: Kenneth Lange @ University of California, Los Angeles 

using StatsBase, Statistics

"""Generates uniform random deviates on the interval (0,1)."""
function uniform_deviate(seed::Vector{Int})
  multiplier = (40014, 40692)
  prime = (2147483563, 2147483399)
  seed[1] = mod(multiplier[1] * seed[1], prime[1])
  seed[2] = mod(multiplier[2] * seed[2], prime[2])
  u = seed ./ prime
  mod(sum(u), 1)
end

n, seed = 100000, [20761807, 58933051]
x = zeros(n)
for i = 1:n
  x[i] = uniform_deviate(seed)
end
println(2 * mean(x), "  ", 12 * var(x)) 
describe(x)
