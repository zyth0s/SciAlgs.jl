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

"""Finds the QR decomposition of A. If A is m x n, then
m <= n should hold."""
function householder_qr(A::Matrix{T}) where T <: Real
  (m, n) = size(A)
  Q = Matrix{T}(I, m, m)
  R = copy(A)
  for j = 1:m
    (c, v) = householder(R[j:m, j])
    R[j:m, :] = R[j:m, :] - (c * v) * (v' * R[j:m, :])
    Q[:, j:m] = Q[:, j:m] - (Q[:, j:m] * v) * (c * v)'
  end
  return(Q, R)
end

(m, n) = (100, 200);
A = randn(m, n);
(Q, R) = householder_qr(A);
println(norm(A - Q * R, 2)," ",norm(Q' * Q - I, 2))
