# Author: Kenneth Lange @ University of California, Los Angeles 

import SciAlgs.TheBook: perfectpower, powermod

"""Tests whether n is prime by the Miller-Rabin method."""
function prime_test(n::T, k::Int) where T <: Integer
  if n == 2
    return prime = true
  end
  if (n < 2) || iseven(n) || perfectpower(n)
    return prime = false
  end
  s = trailing_zeros(n - 1)
  d = div(n - 1, 2^s)
  for trial = 1:k
    a = rand(2:(n - 2))
    x = powermod(a, d, n) # x = a^d mod n
    if x == 1 # n passes
      continue
    end
    t = s
    while x != n - 1 
      t = t - 1
      if t <= 0
        return prime = false
      end
      x = rem(x^2, n) 
      if x == 1 # 1 encountered before -1 so n fails
        return prime = false
      end
    end
  end
  prime = true
end

