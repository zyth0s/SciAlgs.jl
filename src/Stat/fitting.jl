
# Perform least squares minimization 
# min ∑ᵢ [ yᵢ - f(xᵢ,a,b) ]²
# with linear equation f(xᵢ,a,b) = a + bxᵢ
# a: y-intercept
# b: slope
# r: correlation factor ∈ [-1,1]

function lsq_fit(x,y)
   @assert length(x) == length(y)
   N   = length(x) 
   ∑x  = sum(x)
   ∑y  = sum(y)
   x̄ = ∑x / N
   ȳ = ∑y / N
   ∑x² = sum(x.^2)
   ∑y² = sum(y.^2)
   ∑xy = sum( x .* y) # ≡ dot(x,y); but needs LinearAlgebra
   SSxx  = ∑x² - N * x̄^2
   SSyy  = ∑y² - N * ȳ^2
   SSxy  = ∑xy - N * x̄*ȳ
   #covxx =  SSxx / N
   #covyy =  SSyy / N
   #covxy =  SSxy / N

   #b = ( N * ∑xy - ∑x * ∑y ) / ( N * ∑x² - (∑x)^2 )
   #r = ( ∑xy - ∑x * ∑y / N ) / sqrt( (∑x² - (∑x)^2/N) * (∑y² - (∑y)^2/N) )
   a = ( ∑y * ∑x² - ∑x * ∑xy ) / ( N * ∑x² - (∑x)^2 ) 
   b = SSxy/ SSxx
   r = SSxy^2 / (SSxx * SSyy)
   # Standard errors
   s = sqrt( (SSyy - b * SSxy) / (N-2) )
   SEa = s * sqrt( 1/N + x̄^2/SSxx )
   SEb = s / sqrt(SSxx)
   a, b, r, SEa, SEb
end

# Theil's method: 10.1021/ed082p1472.2
function theil_method(x,y)
   N = length(x)
   # remove middle point if N is odd
   if isodd(N)
      deleteat!(x,div(N+1,2))
      deleteat!(y,div(N+1,2))
   end
   n = div(N,2)
   b = [(y[n+i] - y[i])/(x[n+i] - x[i]) for i in 1:n]
   sort!(b)
   b = 0.5 * (b[div(n,2)] + b[div(n+2,2)]) # median
   a = [ y[i] - b*x[i] for i in 1:n]
   sort!(a)
   a = 0.5 * (a[div(n,2)] + a[div(n+2,2)]) # median
   a, b
end

function test_lsq_fit()
   #         xᵢ  yᵢ
   table = [ -3 -15;
              1   1;
              5  17;
              8  29;
             12  45]

   x = table[:,1]
   y = table[:,2]

   a, b, r, SEa, SEb = lsq_fit(x,y)
   println(" y = $a ± $SEa + ($b ± $SEb) x, r = $r")
   @assert isapprox(a, -3)
   @assert isapprox(b,  4)
end

function test_theil_method()
   #         xᵢ  yᵢ
   table = [ -3 -15;
              1   1;
              5  17;
              8  29;
             12  45]

   x = table[:,1]
   y = table[:,2]

   a, b = theil_method(x,y)
   println(" y = $a + ($b) x")
   @assert isapprox(a, -3)
   @assert isapprox(b,  4)
end

test_lsq_fit()
test_theil_method()

