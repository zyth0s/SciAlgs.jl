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

"""Determines whether n is a perfect power."""
function perfectpower(n::Integer)
  m = Int(floor(log(2, n))) # log base 2 of n
  prime_list = eratosthenes(m) 
  for j in prime_list
    k = Int(round(n^(1/j)))
    if isequal(k^j, n)
      return true
    end
  end
  false
end

"""Delivers all primes <= n."""
function eratosthenes(n::Integer)
  isprime = trues(n)
  isprime[1] = false # 1 is composite
  for i = 2:round(Int, sqrt(n))
    if isprime[i]
      for j = i^2:i:n # all multiples of i < i^2 already composite
        isprime[j] = false
      end
    end
  end
  filter(x -> isprime[x], 1:n) # eliminate composite numbers
end

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

using Primes
for i = 1:20
  n = rand(1:10000)
  prime = prime_test(n, 20)
  println(prime," ",isprime(n))
end
