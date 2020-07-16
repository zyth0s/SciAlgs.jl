using LinearAlgebra

"""Implements Karmarker's method of linear programming with the Dinklebach
maneuver.  The program should enter in standard form"""
function karmarkar(A::Matrix{T}, b::Vector{T}, c::Vector{T}, 
  x::Vector{T}) where T <: Real
  A = [A -b] # conversion to homogeneous standard form
  c = [c; zero(T)]
  x = [x; one(T)]
  (m, n) = size(A) 
  Bt = [A; ones(T, n)']' # transposed constraint matrix
  g = zeros(T, n) # gradient in projective space
  cost = dot(c, x)
  for iteration = 1:1000
    for j = 1:n
      Bt[j, 1:m] = A[1:m, j] * x[j]
      g[j] = c[j] * x[j]
    end
    g[n] = g[n] - cost
    y = - Bt \ g # projects -g onto column space of Bt  
    z = - g .- Bt * y # projected steepest descent direction
    z = z / norm(z)
    u = one(T) / n .+ (one(T) / (3 * n)) .* z # new point
    x = u .* x # revert to original coordinates
    x = x ./ x[n]
    newcost = dot(c, x)
    if abs(cost - newcost) < 1e-8 || newcost < -1e6
      return (iteration, x[1:n - 1])
    else
      cost = newcost
    end 
  end
  return (1000, x[1:n - 1])
end

"""Orchestrates Karmarkar's method."""
function karmarkar_program(A::Matrix{T}, b::Vector{T}, c::Vector{T}, 
  tol::T) where T <: Real
  (m, n) = size(A)
  x = rand(n)
  v = b .- A * x
  A1 = [A v]
  c1 = [zeros(T, n); one(T)]
  push!(x, one(T))
  (iterations, x) = karmarkar(A1, b, c1, x) # first phase
  if x[end] > tol
    return ("infeasible", iterations, Inf, x)
  else
    pop!(x)
    (iterations, x) = karmarkar(A, b, c, x) # second phase
    if dot(c, x) < -1e6
      return ("unbounded", iterations, -Inf, x)
    else 
      return ("solvable", iterations, dot(c, x), x)
    end
  end
end

c = [-2.0; -3; -4; 0; 0]; # Wikipedia example
A = [ 3.0 2 1 1 0; 2 5 3 0 1 ];
b = [ 10.0; 15 ];
tol = 1e-5;
(status, iterations, cost, x) = karmarkar_program(A, b, c, tol)
