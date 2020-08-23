# Author: Kenneth Lange @ University of California, Los Angeles 

import SciAlgs.TheBook: eratosthenes

"""Determines whether n is a perfect power."""
function perfectpower(n::Integer)
  m = Int(floor(log(2, n))) # integer part of log base 2 of n
  prime_list = eratosthenes(m) 
  for j in prime_list
    k = Int(round(n^(1 / j)))
    if isequal(k^j, n)
      return true
    end
  end
  false
end
