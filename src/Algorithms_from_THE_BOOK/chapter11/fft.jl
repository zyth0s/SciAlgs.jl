# Author: Kenneth Lange @ University of California, Los Angeles 

using LinearAlgebra

"""Computes the fast Fourier transform of the vector f. ln is the log 
base 2 of the number n of entries of f.  Set invert equal to true for 
the inverse finite Fourier transform and to false for the finite 
Fourier transform.  The  transform is returned in f."""
function FFT(f::Vector{T}, ln::Int, invert::Bool) where T <: Complex
  n = 2^ln
  half_n = div(n, 2)
  j = 1
  for i = 1:(n - 1)
    if i < j
      f[i], f[j] = f[j], f[i]
    end
    k = half_n
    while k < j
      j = j - k
      k = div(k, 2)
    end
    j = j + k
  end
  for l = 1:ln
    k = 2^l
    half_k = div(k, 2)
    u = one(T) 
    w = exp(-pi * im / half_k)
    w = invert ? conj(w) : w
    for j = 1:half_k
      for i = j:k:n
        m = i + half_k
        t = f[m] * u
        f[m] = f[i] - t
        f[i] = f[i] + t
      end
      u = u * w
    end
  end
  invert ? f ./ n : f
end

ln = 4
n = 2^ln
f = randn(n)
f = f .+ 0.0*im
g = copy(f)
f = FFT(f, ln, false)
f = FFT(f, ln, true)
println(norm(f - g))
