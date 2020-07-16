using StatsBase, LinearAlgebra

"""Generates n random deviates according a given mass distribution."""
function discrete_deviate(mass::Vector{T}, n::Int) where T <: Real
  categories = length(mass)
  prob = mass / sum(mass)
  x = zeros(Int, n)
  for trial = 1:n
    u = rand(T)
    s = zero(T)
    for i = 1:categories
      s = s + prob[i]
      if u <= s
        x[trial] = i
        break
      end
    end
  end
  return x
end

(categories, n) = (5, 10000)
mass = rand(categories);
x = discrete_deviate(mass, n);
describe(x)
m = dot([i for i =1:categories], mass / sum(mass)) # match means

