
import Formatting: printfmt
import SpecialFunctions: erf

# Calculates the nodes and weights of the Gauss-Chebyshev 2nd quadrature
# xᵢ = -cos(iθ) where θ = π/(n+1)
# wᵢ 1/√(1-xᵢ²) = θ sin²(iθ) 1/√(1-xᵢ²) = θ sin²(iθ) csc(iθ) = θ sin(iθ) 
function gauss_chebyshev2nd(n)
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  θ::BigFloat = π / (n + 1)
  for i in 1:n
    x[i] = -cos(i*θ)
    w[i] = θ * sin(i*θ)
  end
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


function test_gauss_chebyshev2nd()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Gauss-Chebyshev 2nd kind quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  #for n in [10,20,30,40,50,60,70,80,90]
  for n in 2:10
    x, w = gauss_chebyshev2nd(n)
    println(x)
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
      test_gauss_chebyshev2nd()
    end
  else
      test_gauss_chebyshev2nd()
  end
end
