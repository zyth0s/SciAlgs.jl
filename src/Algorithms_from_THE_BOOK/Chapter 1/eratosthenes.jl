"""Delivers all primes <= n."""
function eratosthenes(n::Integer)
  isprime = trues(n)
  isprime[1] = false
  for i = 2:round(Int, sqrt(n))
    if isprime[i]
      for j = i^2:i:n # all multiples of i < i^2 already composite
        isprime[j] = false
      end
    end
  end
  return filter(x -> isprime[x], 1:n) # eliminate composite numbers
end

prime_list = eratosthenes(100)

