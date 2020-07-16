using LinearAlgebra

"""Constructs a Householder transformation P = I - c u u^t
that maps x into ||x|| e_1."""
function householder(x::Vector{T}) where T <: Real
  a = norm(x[2:end])^2
  u = copy(x)
  u[1] = one(T)
  if a <= zero(T)
    c = zero(T)
  else
    b = sqrt(x[1]^2 + a)
    if x[1] <= zero(T) # avoid roundoff
      u[1] = x[1] - b
    else
      u[1] = - a / (x[1] + b)
    end
    c = 2u[1]^2 / ( a + u[1]^2)
    u = u / u[1]
  end
  return (c, u)
end

"""Reduces the symmetric matrix A to tridiagonal matrix M. 
Forming O' * M * O gives back A."""
function symtridiag(A::Symmetric{T}) where T <: Real
  n = size(A, 2)
  O = Matrix{T}(I, n, n)
  for k = 1:(n - 2)
    (c, u) = householder(A[k + 1:n, k])
    u = [zeros(T, k); u]
    v = A * u
    q = c * v - ((c^2 / 2) * dot(u, v)) * u 
    A = A - q * u' - u * q'
    O = O - (O * u) * (c * u)'
  end
  return (SymTridiagonal(Symmetric(A)), O)
end

n = 100
A = Symmetric(randn(n, n));
B = copy(A);
(T, O) = symtridiag(A);
println(norm(I - O * O'),"  ",norm(B - O * T * O'))
