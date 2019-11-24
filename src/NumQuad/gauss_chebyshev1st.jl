
import Formatting: printfmt
import SpecialFunctions: erf

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

# Tested against https://keisan.casio.com/exec/system/1281438499
function test_gauss_chebyshev1st()
  a = BigFloat(-1); b = BigFloat(1)
  for n in 1:4
    x, w = gauss_chebyshev1st(a,b,n)
    #@assert isapprox(sum(w),b-a) "$(sum(w)) ≠ $(b-a)"
    if n == 1
      @assert isapprox(x[1], (a+b)/2) "$(x[1]) ≠ $((a+b)/2)"
      #@assert isapprox(w[1], ) "$(w[1]) ≠ $(b-a)"
    elseif n == 2
      @assert isapprox(x[1], -0.7071067811865475244008) "$(x[1]) ≠ -0.7071067811865475244008"
      @assert isapprox(x[2],  0.7071067811865475244008) "$(x[2]) ≠  0.7071067811865475244008"
      @assert isapprox(w[1],  1.110720734539591561754 ) "$(w[1]) ≠  1.110720734539591561754"
      @assert isapprox(w[2],  1.110720734539591561754 ) "$(w[1]) ≠  1.110720734539591561754"
    elseif n == 3
      @assert isapprox(x[1], -0.8660254037844386467637) "$(x[1]) ≠ -0.8660254037844386467637"
      @assert isapprox(x[2]+1.0,  1.000000000000000000) "$(x[2]) ≠  0.0000000000000000000000"
      @assert isapprox(x[3],  0.8660254037844386467637) "$(x[2]) ≠  0.8660254037844386467637"
      @assert isapprox(w[1],  0.523598775598298873077) "$(w[1]) ≠ 0.523598775598298873077"
      @assert isapprox(w[2],  1.047197551196597746154) "$(w[1]) ≠ 1.047197551196597746154"
      @assert isapprox(w[3],  0.523598775598298873077) "$(w[1]) ≠ 0.523598775598298873077"
    elseif n == 4
      @assert isapprox(x[1], -0.9238795325112867561282) "$(x[1]) ≠ -0.9238795325112867561282"
      @assert isapprox(x[2], -0.3826834323650897717285) "$(x[2]) ≠ -0.3826834323650897717285"
      @assert isapprox(x[3],  0.3826834323650897717285) "$(w[1]) ≠  0.3826834323650897717285"
      @assert isapprox(x[4],  0.9238795325112867561282) "$(w[1]) ≠  0.9238795325112867561282"
      @assert isapprox(w[1], 0.3005588649421731353571) "$(x[2]) ≠ 0.3005588649421731353571"
      @assert isapprox(w[2], 0.7256132880348577535144) "$(w[1]) ≠ 0.7256132880348577535144"
      @assert isapprox(w[3], 0.7256132880348577535144) "$(w[1]) ≠ 0.7256132880348577535144"
      @assert isapprox(w[4], 0.3005588649421731353571) "$(w[1]) ≠ 0.3005588649421731353571"
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
