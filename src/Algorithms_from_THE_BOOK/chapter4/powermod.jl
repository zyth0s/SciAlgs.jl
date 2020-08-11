# Author: Kenneth Lange @ University of California, Los Angeles 

"""Computes b^m mod n."""
function powermod(b::T, m::T, n::T) where T <: Integer
  if m == zero(T)
    return one(T)
  else
    t = powermod(b, div(m, 2), n)
    t = t^2
    if rem(m, 2) == one(T)
      t = t * b
    end
    return rem(t, n)
  end
end

powermod(2, 6, 5)

