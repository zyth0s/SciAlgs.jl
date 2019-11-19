
import Formatting: printfmt
import SpecialFunctions: erf

#using QuadGK

# This routine computes the points and weights of the Gauss-Chebyshev
# quadrature of the second kind using the Perez-Jorda transformation.
# "A simple, efficient and more reliable scheme for automatic numerical
# integration" CPC eqs. 20-21.
function perez_jorda(n)
  # Limited to [-1,1]
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  for i in 1:n
    dn = n+1
    θ = i * π / dn
    st = sin(θ)
    ct = cos(θ)
    w[i] = 16.0 / 3.0 / dn * st * st * st * st
    x[i] = -1.0 + 2.0 * i / dn - 2.0 / π * (1 + 2.0/3.0 * st * st) * ct * st
  end
  if n == 1
    w[1] = 2.00
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

function test_perez_jorda()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Perez-Jorda quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = perez_jorda(n)
    numint = sum( map(f, x) .* w)
    error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
    bar = "▋"^(  digits + ceil(Int,get_log10(abs(error),digits))     )
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}," * 
             " error is {:$(digits+4).$(digits)f} " *
             bar * "\n",
            n,numint, error)
  end
  #numint, error = quadgk(f,a,b,rtol=exp10(-digits))
  #printfmt("    Integral with QuadGK is    {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",numint, error)
  #error = numint - 0.5*√(π)*(erf(b) - erf(a) )
  #printfmt("    Error of QuadGK against ref    {:$(digits+4).$(digits)f}\n", error)
end

if ! isinteractive()
  if abspath(PROGRAM_FILE) == @__FILE__
    setprecision(100) do         # standard double precision is 53
      test_perez_jorda()
    end
  else
      test_perez_jorda()
  end
end
