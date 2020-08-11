# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra, SparseArrays

""" Solves the equation A*x = b by the conjugate gradient method. 
The matrix A is assumed positive definite; tol is tolerance for 
testing convergence."""
function conjugategradients(A::AbstractMatrix{T}, b::Vector{T}, 
  tol::T) where T <: Real
#
  n = size(A, 1)
  if n == 1
    return b / A[1, 1]
  else
    x = zeros(T, n) # solution vector starts at 0
    v = copy(b) # search direction
    r = copy(b) # residual b - A * x
    s = norm(r)^2 # residual sum of squares
    for j = 1:n # perform the conjugate gradient steps 
      av = A * v
      c = dot(v, av)
      if c <= zero(T)
        return x
      end
      t = s / c # step length
      x = x + t * v
      r = r - t * av
      d = norm(r)^2
      if sqrt(d) < tol # convergence test
        return x
      end
      v = r + (d / s) * v
      s = d
    end
    return x
  end
end

n, tol = 100, 1e-5
A = sprandn(n, n, .01)
rho = opnorm(A, 1)
A = 0.5 * (A + A') + rho * I
x = randn(n)
b = A * x
y = conjugategradients(A, b, tol)
norm(x - y)
