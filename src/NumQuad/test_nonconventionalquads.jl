

using Formatting: printfmt
using LinearAlgebra
using SpecialFunctions
#using QuadGK

include("maxwell.jl")
include("multiexp.jl")

digits = trunc(Int,log10( 2.0^precision(Float64)))

function maxwell_w(r,kind=1)
  @assert 0 <= kind
  1/ ( r^kind * exp(-r^2) )
end

function maxwell_tester(f,analytic_result)

  for kind in 0:2
    println("Testing with Maxwell polynomials of kind $kind")
    for n in 10:10:90
      r, w = maxwellpts(n,kind)
      numint = sum( @. w * f(r) * maxwell_w(r,kind))
      error = numint - analytic_result
      printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n, numint, error)
    end
    println("")
  end
end

function multiexp_tester(f,analytic_result)

  println("Testing with Multiexp grid")
  for n in 10:10:90
    r, w = multiexp(n)
    numint = sum( @. w * f(r) / log(r)^2 )
    error = numint - analytic_result
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n, numint, error)
  end
  println("")
end

function fejer_tester(f,analytic_result)

  println("Testing with composite Fejer quadrature rule")
  for n in 1:1:10
    pwmd = multidomain_quadrature(n,n,-1,1)
    r = pwmd[:,1]
    w = pwmd[:,2]
    numint = sum( @. w * f(r) )
    error = numint - analytic_result
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n*n, numint, error)
  end
  println("")
end

function test_maxwell()
  a = 0; b = 30
  println("Integrating f(r) = r*exp(-r)")
  println("================================")
  f(r) = r*exp(-r)
  analytic_result = ( -exp(-b)*(1+b) + exp(-a)*(1+a) )
  maxwell_tester(f,analytic_result)

  println("Integrating g(x) = exp(-x^2)")
  println("================================")
  g(x) = exp(-(x)^2)
  analytic_result = 0.5*√(π)*(erf(b) - erf(a) )
  maxwell_tester(g,analytic_result)
end

function test_multiexp()
  a = 0; b = 1
  println("Integrating f(r) = r*exp(-r)")
  println("================================")
  f(r) = r*exp(-r)
  analytic_result = ( -exp(-b)*(1+b) + exp(-a)*(1+a) )
  multiexp_tester(f,analytic_result)

  println("Integrating g(x) = exp(-x^2)")
  println("================================")
  g(x) = exp(-(x)^2)
  analytic_result = 0.5*√(π)*(erf(b) - erf(a) )
  multiexp_tester(g,analytic_result)
end

function test_fejer()
  a = -1; b = 1
  println("Integrating f(r) = r*exp(-r)")
  println("================================")
  f(r) = r*exp(-r)
  analytic_result = ( -exp(-b)*(1+b) + exp(-a)*(1+a) )
  fejer_tester(f,analytic_result)

  println("Integrating g(x) = exp(-(x+1)^2)")
  println("================================")
  g(x) = exp(-(x+1)^2)
  analytic_result = 0.5*√(π)*(erf(1+b) - erf(1+a) )
  fejer_tester(g,analytic_result)
end


test_maxwell()
test_multiexp()
test_fejer()
