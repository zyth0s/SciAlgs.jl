
import Formatting: printfmt
import SpecialFunctions: erf

# Adapted algorithm of Homer Reid to [a,b]
# Learn more from "Fast Construction of the Fejer and Clenshaw-Curtis Quadrature"
# by Joerg Waldvogel eqs. (2.4-2.5)
function clenshaw_curtis(a::BigFloat,b::BigFloat,n)
  if n < 1
    error("At least one quadrature point should be requested")
  end
  #x = zeros(n)
  x1 = (b - a)*0.5
  x2 = (b + a)*0.5
  nn = n-1
  nnh = nn/2
  x = @. -cos((0:nn)*π/nn) * x1 + x2
  w = ones(BigFloat,n)
  for i in 0:nn
    θ = i * π / nn
    for j = 1:nnh
      if (2*j) == nn
        bj = 1.0
      else
        bj = 2.0
      end
      w[i+1] -= bj * cos(2.0 * j * θ) / (4*j*j - 1)
    end
    #w[i+1] *= (b-a) / nn
    w[i+1] *= (b-a) / 2.0 / nn
  end
  #w[1] /= nn
  w[2:nn] *= 2.0
  #w[n] /= nn
  #w[1] /= nn
  #w[2:n-1] *= 2.0 / nn
  #w[n] /= nn
  if n < 2 # == 1
    w[1] = b-a
    x[1] = x2
  end
  #for i in eachindex(x)
  #  println(x[i], "  ",w[i])
  #end
  x, w
end

function get_log10(x,digits)
  l = log10(x)
  if isinf(l)
    return -digits
  else
    return l
  end
end

function test_clenshaw_curtis2()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Clenshaw-Curtis quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = clenshaw_curtis(a,b,n)
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
      #test_clenshaw_curtis()
      test_clenshaw_curtis2()
    end
  end
end
