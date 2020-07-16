using LinearAlgebra

"""Performs matrix multiplication by random outer products."""
function randommultiply(A::AbstractMatrix{T}, B::AbstractMatrix{T}, 
  r::Int) where T <: Real
  u = sqrt.(sum(abs2, A, dims = 1)) # norms of columns of A
  v = sqrt.(sum(abs2, B, dims = 2)) # norms of rows of B
  p = vec(u) .* vec(v)
  p = p / sum(p) 
  s = discrete_deviate(p, r) # random sample from p
  AA = zeros(T, size(A, 1), 0)
  BT = zeros(T, size(B, 2), 0)
  for i = 1:r
    AA = [AA (p[s[i]]^(-1) .* A[:, s[i]])]
    BT = [BT B[s[i], :]]
  end  
  return (AA * BT') ./ r # instead of summing outer products
end
  
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

(m, n, p, r) = (1000, 100, 1000, 100);
A = randn(m, n);
B = randn(n, p);
@time C = A * B;
@time D = randommultiply(A, B, r);
println(norm(C - D) / (m * p))
