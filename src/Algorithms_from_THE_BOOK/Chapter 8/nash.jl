using LinearAlgebra

"""Constructs a rotation that zeros out the offdiagonal 
entry b of a 2 x 2 symmetric matrix."""
function jacobi_rotation(a::T, b::T, c::T) where T <: Real
  ratio = (c - a) / (2b)
  if ratio >= zero(T)
    tangent = one(T) / (ratio + sqrt(one(T) + ratio^2))
  else
    tangent = -one(T) / (abs(ratio) + sqrt(one(T) + ratio^2))
  end
  cosine = one(T) / sqrt(one(T) + tangent^2)
  sine = tangent * cosine
  return (cosine, sine)
end

"""Extracts the singular value decomposition U * diag(d) * V' of a 
matrix M by the method of Nash. M is destroyed in the process.
The number of rows of M should exceed the number of columns."""
function nash(M::Matrix{T}, tol::T) where T <: Real
  (m, n) = size(M)
  Mnorm = norm(M)
  V = Matrix{T}(I, n, n)
  converge = Mnorm + one(T)
  while abs(converge) > tol * (Mnorm + one(T)) # perform rotations
    converge = zero(T)
    for i = 1:n
      for j = (i + 1):n
        (u, v) = (M[:, i], M[:, j])
        (a, b, c) = (norm(u)^2, dot(u, v), norm(v)^2)
        converge = max(converge, abs(b))
        (cosine, sine) = jacobi_rotation(a, b, c)
        M[:, i] = cosine * u - sine * v
        M[:, j] = sine * u + cosine * v
        (u, v) = (V[:, i], V[:, j])
        V[:, i] = cosine * u - sine * v
        V[:, j] = sine * u + cosine * v
      end
    end
  end
  d = zeros(T, n) # find singular values and left vectors
  for i = 1:n
    d[i] = norm(M[:, i])
    M[:, i] = M[:, i] ./ d[i] # column scaling
  end
  return (M, d, V) # singular values not in decreasing order
end

(m, n, tol) = (100, 50, 1e-12);
M = randn(m, n);
Msave = copy(M);
@time (U, d, V) = nash(M, tol);
norm(Msave - U * (Diagonal(d) * V'))




