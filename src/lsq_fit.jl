
# Perform least squares minimization 
# min [ y - yᵢ ]²
# with linear equation y = mx + b
# m: slope
# b: y-intercept

function lsq_fit(x,y)
   N = length(x) 

   ∑x  = sum(x)
   ∑y  = sum(y)
   ∑x² = sum(x.^2)
   ∑y² = sum(y.^2)
   ∑xy = sum( x .* y) # ≡ dot(x,y); but needs LinearAlgebra

   m = ( N * ∑xy - ∑x * ∑y ) / ( N * ∑x² - (∑x)^2 )
   b = ( ∑y * ∑x² - ∑x * ∑xy ) / ( N * ∑x² - (∑x)^2 ) 

   r = ( ∑xy - ∑x * ∑y / N ) / sqrt( (∑x² - (∑x)^2/N) * (∑y² - (∑y)^2/N) )
   m, b, r
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

   m, b, r = lsq_fit(x,y)
   @assert m ==  4
   @assert b == -3
   println(" y = $m x + ($b), r = $r")
end

test_lsq_fit()

