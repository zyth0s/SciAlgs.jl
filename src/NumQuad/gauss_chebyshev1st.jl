
import Formatting: printfmt
import SpecialFunctions: erf

import SciAlgs.NumQuad: get_log10

# Calculates the nodes and weights of the Gauss-Chebyshev 1st quadrature
# xᵢ = -cos(θᵢ) where θᵢ = (2i-1)*π/2/n
# wᵢ √(1-xᵢ²) = π/n √(1-xᵢ²) = π/n sin(θᵢ) 
# scaled to [a,b]
function gauss_chebyshev1st(a::BigFloat,b::BigFloat,n)
  x1 = (b-a)*0.5
  x2 = (b+a)*0.5
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  for i in 1:n
    θ = (i-0.5)*π/n
    x[i] = -cos(θ)*x1 + x2
    w[i] = π/n * x1 * sin(θ)
  end
  if n == 1
    x[1] = x2
    w[1] = b-a
  end
  #for i in eachindex(x)
  #  println(x[i], "  ",w[i])
  #end
  x, w
end

function test_gauss_chebyshev1st2()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Gauss-Chebyshev 1st kind quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = gauss_chebyshev1st(a,b,n)
    numint = sum( map(f, x) .* w)
    error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
    bar = "▋"^(  digits + ceil(Int,get_log10(abs(error),digits))     )
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}," * 
             " error is {:$(digits+4).$(digits)f} " *
             bar * "\n",
            n,numint, error)
  end
  println("------------------------------------------------------")
  println("Gauss-Chebyshev 1st kind quadrature")
  println("Testing integration of f(r) = 1/√(1 - r^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  g(x) = 1/ √(1 - x^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = gauss_chebyshev1st(a,b,n)
    numint = sum( map(g, x) .* w)
    error = numint - π
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
      test_gauss_chebyshev1st()
      test_gauss_chebyshev1st2()
    end
  end
end
