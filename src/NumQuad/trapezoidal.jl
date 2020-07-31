

import Formatting: printfmt
import SpecialFunctions: erf

import SciAlgs.NumQuad: get_log10

# Calculates the nodes and weights of the trapezoidal quadrature
# x_i = 
# w_i =
function trapezoidal(n)
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  u = BigFloat(-1)
  Δu =  BigFloat(2) / (n - 1)
  w = fill(Δu,n)
  w[1]   = Δu / 2
  w[end] = Δu / 2
  for i in 1:n
    x[i] = u
    u += Δu
  end
  x, w
end

function test_trapezoidal()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Trapezoidal quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = trapezoidal(n)
    numint = sum( map(f, x) .* w)
    error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
    bar = "▋"^(  digits + ceil(Int,get_log10(abs(error),digits))     )
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}," * 
             " error is {:$(digits+4).$(digits)f} " *
             bar * "\n",
            n,numint, error)
  end
end

if ! isinteractive()
  if abspath(PROGRAM_FILE) == @__FILE__
    setprecision(100) do         # standard double precision is 53
      test_trapezoidal()
    end
  end
end
