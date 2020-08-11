# Author: Kenneth Lange @ University of California, Los Angeles 

""" Computes the greatest common divisor of c and d. """
function euclid(m::T, n::T) where T <: Integer
  (a, b) = (m, n)
  while b != zero(T)
    a, b = b, rem(a, b)
  end
  a
end

gcd = euclid(600, 220)
