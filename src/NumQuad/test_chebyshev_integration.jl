
using Formatting: printfmt
using LinearAlgebra
using SpecialFunctions
#using QuadGK
using FastGaussQuadrature

digits = trunc(Int,log10( 2.0^precision(Float64)))

function cheb_w(r,kind=1)
  @assert 0 < kind < 5
  if kind == 1
    return √(1-r^2)
  elseif kind == 2
    return 1/√(1-r^2)
  elseif kind == 3
    return √((1-r)/(1+r))
  else kind == 4
    return √((1+r)/(1-r))
  end
end

#function example1()
#
#  f(r) = r*exp(-r)
#
#  a = -1.0 # gausschebyshev is fixed to [-1,1]
#  b =  1.0 # " " "
#
#  for kind in 1:4
#    println("Testing with Chebyshev polynomials of kind $kind")
#    for n in 10:10:90
#      r, w = gausschebyshev(n,kind)
#      numint = sum( @. w * f(r) * cheb_w(r,kind))
#      error = numint - ( -exp(-b)*(1+b) + exp(-a)*(1+a) )
#      printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n, numint, error)
#    end
#    println("")
#  end
#end

function cheb_tester(f,analytic_result)

  for kind in 1:4
    println("Testing with Chebyshev polynomials of kind $kind")
    for n in 10:10:90
      r, w = gausschebyshev(n,kind)
      numint = sum( @. w * f(r) * cheb_w(r,kind))
      error = numint - analytic_result
      printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n, numint, error)
    end
    println("")
  end
end

function test()
  a = -1; b = 1 # Remember that Chebyshev polynomials are defined in the range [-1,1]
  #example1()
  # Reproduce the test
  # https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/48931/versions/4/previews/INT_TLBX/html/GCQuadVerify.html
  println("Integrating f(r) = r*exp(-r)")
  println("================================")
  f(r) = r*exp(-r)
  analytic_result = ( -exp(-b)*(1+b) + exp(-a)*(1+a) )
  cheb_tester(f,analytic_result)
  # Observation: Gauss-Chebyshev 1st/2nd kind has a weighting that has quite
  # the opposite shape of the integrand. That is why the integration is so
  # inaccurate.

  println("Integrating g(x) = exp(-(x+1)^2)")
  println("================================")
  g(x) = exp(-(x+1)^2)
  analytic_result = 0.5*√(π)*(erf(1+b) - erf(1+a) )
  cheb_tester(g,analytic_result)
end

#test()
