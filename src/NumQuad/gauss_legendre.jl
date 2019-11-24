
import Formatting: printfmt
import SpecialFunctions: erf


function gauss_legendre(x1::BigFloat, x2::BigFloat, n)
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  eps = 3e-16
  m = div(n+1, 2)
  xm = 0.5 * (x2 + x1)
  xl = 0.5 * (x2 - x1)
  for i in 1:m
    z = cos(π * (i - 0.25) / (n + 0.5))
    @label label1
    p1 = 1.0
    p2 = 0.0
    for j in 1:n
      p3 = p2
      p2 = p1
      p1 = ((2.0*j - 1.0) * z * p2 - (j - 1.0)*p3) / j
    end
    pp = n * (z * p1 - p2) / (z * z - 1.0)
    z1 = z
    z = z1 - p1 / pp
    if abs(z - z1) > eps
      @goto label1
    end
    x[i] = xm - xl * z
    x[n+1-i] = xm + xl * z
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp)
    w[n+1-i] = w[i]
  end
  #for i in eachindex(x)
  #  println(x[i], "  ",w[i])
  #end
  x, w
end

function test_gauss_legendre()
  a = BigFloat(-1); b = BigFloat(1)
  for n in 1:6
    x, w = gauss_legendre(a,b,n)
    @assert isapprox(sum(w),b-a) "$(sum(w)) ≠ $(b-a)"
    eps = 1e-5
    if n == 1
      @assert isapprox(x[1], (a+b)/2) "$(x[1]) ≠ $((a+b)/2)"
      #
      @assert isapprox(w[1], b-a) "$(w[1]) ≠ $(b-a)"
    elseif n == 2
      @assert isapprox(x[1], -0.577350,atol=eps) "$(x[1]) ≠ $(-0.577350)"
      @assert isapprox(x[2], +0.577350,atol=eps) "$(x[1]) ≠ $(+0.577350)"
      #
      @assert isapprox(w[1], 1.000000,atol=eps) "$(w[1]) ≠ $(1.000000)"
      @assert isapprox(w[2], 1.000000,atol=eps) "$(w[1]) ≠ $(1.000000)"
    elseif n == 3
      @assert isapprox(x[1], -0.774597,atol=eps) "$(x[1]) ≠ $(-0.774597)"
      @assert isapprox(x[2],  0.000000,atol=eps) "$(x[1]) ≠ $( 0.000000)"
      @assert isapprox(x[3], +0.774597,atol=eps) "$(x[1]) ≠ $(+0.774597)"
      #
      @assert isapprox(w[1], 0.555556,atol=eps) "$(w[1]) ≠ $(0.555556)"
      @assert isapprox(w[2], 0.888889,atol=eps) "$(w[1]) ≠ $(0.888889)"
      @assert isapprox(w[3], 0.555556,atol=eps) "$(w[1]) ≠ $(0.555556)"
    elseif n == 4
      @assert isapprox(x[1], -0.861136,atol=eps) "$(x[1]) ≠ $(-0.861136)"
      @assert isapprox(x[2], -0.339981,atol=eps) "$(x[1]) ≠ $(-0.339981)"
      @assert isapprox(x[3], +0.339981,atol=eps) "$(x[1]) ≠ $(+0.339981)"
      @assert isapprox(x[4], +0.861136,atol=eps) "$(x[1]) ≠ $(+0.861136)"
      #
      @assert isapprox(w[1], 0.347855,atol=eps) "$(w[1]) ≠ $(0.347855)"
      @assert isapprox(w[2], 0.652145,atol=eps) "$(w[1]) ≠ $(0.652145)"
      @assert isapprox(w[3], 0.652145,atol=eps) "$(w[1]) ≠ $(0.652145)"
      @assert isapprox(w[4], 0.347855,atol=eps) "$(w[1]) ≠ $(0.347855)"
    elseif n == 5                                                    
      @assert isapprox(x[1],-0.906180,atol=eps) "$(x[1]) ≠ $(-0.906180)"
      @assert isapprox(x[2],-0.538469,atol=eps) "$(x[1]) ≠ $(-0.538469)"
      @assert isapprox(x[3], 0.000000,atol=eps) "$(x[1]) ≠ $( 0.000000)"
      @assert isapprox(x[4],+0.538469,atol=eps) "$(x[1]) ≠ $(+0.538469)"
      @assert isapprox(x[5],+0.906180,atol=eps) "$(x[1]) ≠ $(+0.906180)"
      #
      @assert isapprox(w[1],0.236927,atol=eps) "$(w[1]) ≠ $(0.236927)"
      @assert isapprox(w[2],0.478629,atol=eps) "$(w[1]) ≠ $(0.478629)"
      @assert isapprox(w[3],0.568889,atol=eps) "$(w[1]) ≠ $(0.568889)"
      @assert isapprox(w[4],0.478629,atol=eps) "$(w[1]) ≠ $(0.478629)"
      @assert isapprox(w[5],0.236927,atol=eps) "$(w[1]) ≠ $(0.236927)"
    elseif n == 6
      @assert isapprox(x[1],-0.932470,atol=eps) "$(x[1]) ≠ $(-0.932470)"
      @assert isapprox(x[2],-0.661209,atol=eps) "$(x[1]) ≠ $(-0.661209)"
      @assert isapprox(x[3],-0.238619,atol=eps) "$(x[1]) ≠ $(-0.238619)"
      @assert isapprox(x[4],+0.238619,atol=eps) "$(x[1]) ≠ $(+0.238619)"
      @assert isapprox(x[5],+0.661209,atol=eps) "$(x[1]) ≠ $(+0.661209)"
      @assert isapprox(x[6],+0.932470,atol=eps) "$(x[1]) ≠ $(+0.932470)"
      #
      @assert isapprox(w[1],0.1713244924,atol=eps) "$(w[1]) ≠ $(0.1713244924)"
      @assert isapprox(w[2],0.360761753 ,atol=eps) "$(w[1]) ≠ $(0.360761753 )"
      @assert isapprox(w[3],0.467914    ,atol=eps) "$(w[1]) ≠ $(0.467914    )"
      @assert isapprox(w[4],0.467914    ,atol=eps) "$(w[1]) ≠ $(0.467914    )"
      @assert isapprox(w[5],0.360761753 ,atol=eps) "$(w[1]) ≠ $(0.360761753 )"
      @assert isapprox(w[6],0.1713244924,atol=eps) "$(w[1]) ≠ $(0.1713244924)"
    end
  end
end

function get_log10(x,digits)
  l = log10(x)
  if isinf(l)
    return -digits
  else
    return l
  end
end

function test_gauss_legendre2()
  a = BigFloat(-1); b = BigFloat(1)
  println("------------------------------------------------------")
  println("Gauss-Legendre quadrature")
  println("Testing integration of f(r) = exp(-(r+1)^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  f(x) = exp(-(x+1)^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = gauss_legendre(a,b,n)
    numint = sum( map(f, x) .* w)
    error = numint - 0.5*√(π)*(erf(1+b) - erf(1+a) )
    bar = "▋"^(  digits + ceil(Int,get_log10(abs(error),digits))     )
    printfmt("    Integral with {:7d} points is {:$(digits+4).$(digits)f}," * 
             " error is {:$(digits+4).$(digits)f} " *
             bar * "\n",
            n,numint, error)
  end
  println("------------------------------------------------------")
  println("Gauss-Legendre quadrature")
  println("Testing integration of f(r) = 1/√(1 - r^2) from $a to $b")
  println("------------------------------------------------------")
  println("NOTE: Testing with $(precision(BigFloat)) bits of precision")
  digits = trunc(Int,log10( 2.0^precision(BigFloat)))
  g(x) = 1/ √(1 - x^2)
  for n in [10,20,30,40,50,60,70,80,90]
    x, w = gauss_legendre(a,b,n)
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
      #test_gauss_legendre()
      test_gauss_legendre2()
    end
  end
end

