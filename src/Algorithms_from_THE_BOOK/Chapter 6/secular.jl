""" Computes the roots of the quadratic a*x^2+b*x+c when |a|>0."""
function quadratic(a::T, b::T, c::T) where T <: Real
  d = b^2 - 4a * c # discriminant
  if d > zero(T) 
    if b >= zero(T)
      r1 = (-b - sqrt(d)) / (2a)
    else
      r1 = (-b + sqrt(d)) / (2a)
    end
    r2 = c / (r1 * a)
    return (r1, r2)
  else
    return (-b + sqrt(d + 0im)) / (2a), (-b - sqrt(d + 0im)) / (2a)
  end
end

"""Finds the root of the surrogate for the secular function."""
function approx_secular(d::Vector{T}, w::Vector{T}, mu::T,
  lambda::T, j::Int) where T <: Real
  n = length(d)
  (f1, f2, c) = (zeros(T, 2), zeros(T, 2), zeros(T, 3))
  for i = 1:j # values and derivatives of the secular sums
    f1[1] = f1[1] + w[i] / (d[i] - lambda)
    f1[2] = f1[2] + w[i] / (d[i] - lambda)^2
  end
  for i = (j + 1):n
    f2[1] = f2[1] + w[i] / (d[i] - lambda)
    f2[2] = f2[2] + w[i] / (d[i] - lambda)^2
  end
  c[3] = mu # coefficients of the surrogate
  if j == 0
    c[2] = f2[2] * (d[j + 1] - lambda)^2
    c[3] = c[3] + f2[1] - f2[2] * (d[j + 1] - lambda)
    return ((c[2] + c[3] * d[j + 1]) / c[3], mu + f1[1] + f2[1])
  elseif j == n
    c[1] = f1[2] * (d[j] - lambda)^2
    c[3] = c[3] + f1[1] - f1[2] * (d[j] - lambda)
    return ((c[1] + c[3] * d[j]) / c[3], mu + f1[1] + f2[1])
  else
    c[1] = f1[2] * (d[j] - lambda)^2
    c[3] = c[3] + f1[1] - f1[2] * (d[j] - lambda)
    c[2] = f2[2] * (d[j + 1] - lambda)^2
    c[3] = c[3] + f2[1] - f2[2] * (d[j + 1] - lambda)
    p = c[3] # coefficients of the quadratic
    q = - (c[1] + c[2] + c[3] * (d[j] + d[j + 1]))
    r = c[1] * d[j + 1] + c[2] * d[j] + c[3] * d[j] *  d[j + 1]
    (root1, root2) = quadratic(p, q, r) # solve the quadratic
    if root1 > d[j] && root1 < d[j + 1]
      return (root1, mu + f1[1] + f2[1])
    else
      return (root2, mu + f1[1] + f2[1])
    end
  end
end

"""Solves the secular equation by the LAPACK method."""
function solve_secular(d::Vector{T}, w::Vector{T}, mu::T, 
  tol::T) where T <: Real
  n = length(d)
  (lambda, value) = (zeros(T, n), zeros(T, n))
  for i = 1:n
    if mu < zero(T) # choose the interval for ith root 
      j = i - 1
    else
      j = i
    end
    if j == 0 # initialize ith root
      lambda[i] = d[1] - one(T) # to the left
    elseif j == n
      lambda[i] = d[n] + one(T) # to the right
    else
      lambda[i] = (d[j] + d[j + 1]) / 2 # midpoint
    end
    for k = 1:10 # iterate until convergence
      (lambda[i], value[i]) = approx_secular(d, w, mu, lambda[i], j)
      if abs(value[i]) < tol
        break
      end
    end
  end
  return (lambda, value) # roots and secular function values
end

(n, mu, tol) = (100, randn(), 1e-10);
w = rand(n);
d = sort(randn(n));
(lambda, value) = solve_secular(d, w, mu, tol);
for i = 1:length(d)
  println(i,"  ",lambda[i],"  ",value[i])
end
