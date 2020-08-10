
using SciAlgs.NumQuad: gauss_chebyshev2nd

import SpecialFunctions: erf

@testset "NumQuad: Gauss-Chebyshev 2nd order" begin

   # (Crude) testing integration of f(r) = exp(-(r+1)^2) from a to b

   a = BigFloat(-1); b = BigFloat(1)
   f(x) = exp(-(x+1)^2)

   for n in [10,20,30,40,50,60,70,80,90]

      x, w = gauss_chebyshev2nd(n)
      numint = sum( map(f, x) .* w)
      error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
      @test error < 1e-3
   end
end


#function test_gauss_chebyshev2nd()
#  a = BigFloat(-1); b = BigFloat(1)
#  println("------------------------------------------------------")
#  println("Gauss-Chebyshev 2nd kind quadrature")
#  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
#  println("------------------------------------------------------")
#  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
#  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
#  f(x) = exp(-(x+1)^2)
#  #for n in [10,20,30,40,50,60,70,80,90]
#  for n in 2:10
#    x, w = gauss_chebyshev2nd(n)
#    println(x)
#    numint = sum( map(f, x) .* w)
#    error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
#    bar = "▋"^(  digits + ceil(Int,get_log10(abs(error),digits))     )
#    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}," * 
#             " error is {:$(digits+4).$(digits)f} " *
#             bar * "\n",
#            n,numint, error)
#  end
#end
