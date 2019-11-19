
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

function test_clenshaw_curtis()
  for i in [1,2,3,4,5,7]
    for a in BigFloat[-4,-3,-2,-1, 0], b in BigFloat[1,2,3,4,5]
      x, w = clenshaw_curtis(a,b,i)
      @assert isapprox(sum(w),b-a) "$(sum(w)) ≠ $(b-a)"
      if i == 1
        @assert isapprox(w[1], b-a) "$(w[1]) ≠ $(b-a)"
      elseif i == 2
        @assert isapprox(w[1], (b-a)/2) "$(w[1]) ≠ $((b-a)/2)"
        @assert isapprox(w[2], (b-a)/2) "$(w[1]) ≠ $((b-a)/2)"
      elseif i == 3
        @assert isapprox(w[1],   (b-a)/6) "$(w[1]) ≠ $(  (b-a)/6)"
        @assert isapprox(w[2], 4*(b-a)/6) "$(w[1]) ≠ $(4*(b-a)/6)"
        @assert isapprox(w[3],   (b-a)/6) "$(w[1]) ≠ $(  (b-a)/6)"
      elseif i == 4
        @assert isapprox(w[1],   (b-a)/18) "$(w[1]) ≠ $(   (b-a)/18)"
        @assert isapprox(w[2], 8*(b-a)/18) "$(w[1]) ≠ $( 8*(b-a)/18)"
        @assert isapprox(w[3], 8*(b-a)/18) "$(w[1]) ≠ $( 8*(b-a)/18)"
        @assert isapprox(w[4],   (b-a)/18) "$(w[1]) ≠ $(   (b-a)/18)"
      elseif i == 5                                                    
        @assert isapprox(w[1],   (b-a)/30) "$(w[1]) ≠ $(   (b-a)/30)"
        @assert isapprox(w[2], 8*(b-a)/30) "$(w[1]) ≠ $( 8*(b-a)/30)"
        @assert isapprox(w[3],12*(b-a)/30) "$(w[1]) ≠ $(12*(b-a)/30)"
        @assert isapprox(w[4], 8*(b-a)/30) "$(w[1]) ≠ $( 8*(b-a)/30)"
        @assert isapprox(w[5],   (b-a)/30) "$(w[1]) ≠ $(   (b-a)/30)"
      elseif i == 7
        @assert isapprox(w[1],  9*(b-a)/630) "$(w[1]) ≠ $(  9*(b-a)/630)"
        @assert isapprox(w[2], 80*(b-a)/630) "$(w[1]) ≠ $( 80*(b-a)/630)"
        @assert isapprox(w[3],144*(b-a)/630) "$(w[1]) ≠ $(144*(b-a)/630)"
        @assert isapprox(w[4],164*(b-a)/630) "$(w[1]) ≠ $(164*(b-a)/630)"
        @assert isapprox(w[5],144*(b-a)/630) "$(w[1]) ≠ $(144*(b-a)/630)"
        @assert isapprox(w[6], 80*(b-a)/630) "$(w[1]) ≠ $( 80*(b-a)/630)"
        @assert isapprox(w[7],  9*(b-a)/630) "$(w[1]) ≠ $(  9*(b-a)/630)"
      end
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
  else
      #test_clenshaw_curtis()
      test_clenshaw_curtis2()
  end
end
