
function Gauss_Easter_algorithm(A)
   # Donald Teets (2019) Gauss’s Computation of the Easter Date, Mathematics Magazine, 92:2, 91-98
   A > 1582 || error("Must be a year in the Gregorian Calendar")
   a = mod(A,19)
   b = mod(A,4)
   c = mod(A,7)
   k = div(A,100)
   p = div(13+8k,25)
   q = div(k,4)
   M = mod(15-p+k-q,30)
   N = mod(4+k-q,7)
   d = mod(19a+M,30)
   e = mod(2b+4c+6d+N,7)

   domPascua = 0
   if d+e < 10
      domPascua = d+e+22
      println("Easter Sunday is the $domPascua th of March.")
   else
      domPascua = d+e-9
      if domPascua == 26
         domPascua = 19
      elseif domPascua == 25 && d == 28 && e == 6
         domPascua = 18
      end
      println("Easter Sunday is the $domPascua th of April.")
   end
end

# O'Beirne, T. H. "How ten divisions lead to Easter". New Scientist, 9 (228): 828. (30 March 1961).
function Meesus_Jones_Butcher_Easter_algorithm(x)
   x > 1582 || error("Must be a year in the Gregorian Calendar")
   a = mod(x,19)
   b, c = divrem(x,100)
   d, e = divrem(b,4)
   g = div(8b+13,25)
   h = mod(19a+b-d-g+15,30)
   i, k = divrem(c,4)
   l = mod(2e +2i-h-k+32,7) # warn: 2e+2 is interpreted as 200
   m = div(a+11h+19l,433)
   n = div(h+l+7m+90,25)
   p = mod(h+l-7m+33n+19,32)
   println("In the $x th year AD of the Gregorian calendar,") 
   println("Easter Sunday is the $p th day of the $n th month.")
end

Gauss_Easter_algorithm(2020)
println("")
Meesus_Jones_Butcher_Easter_algorithm(2020)
