

using Formatting: printfmt
using LinearAlgebra
#using QuadGK
using FastGaussQuadrature

#include("gauss_laguerre.jl")

digits = trunc(Int,log10( 2.0^precision(Float64)))


f(r) = exp(-r^2) 

a = 0.0 # gausslaguerre is fixed to [0,Inf]
b = Inf # " " "

function test_laguerre()
   for n in 10:10:90
     r, w = gausslaguerre(n,0.0)
     #numint = dot( w, f.(r))
     numint = sum( @. w * f(r) * exp(r))
     error = numint - 0.5*√π # 1 # 3628800.0 # 1 
     printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}, error is {:$(digits+4).$(digits)f}\n",n, numint, error)
   end
end
