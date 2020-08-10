
import SciAlgs.NumQuad: get_log10

# Calculates the nodes and weights of the Gauss-Chebyshev 2nd quadrature
# xᵢ = -cos(iθ) where θ = π/(n+1)
# wᵢ 1/√(1-xᵢ²) = θ sin²(iθ) 1/√(1-xᵢ²) = θ sin²(iθ) csc(iθ) = θ sin(iθ) 
function gauss_chebyshev2nd(n)
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  θ::BigFloat = π / (n + 1)
  for i in 1:n
    x[i] = -cos(i*θ)
    w[i] = θ * sin(i*θ)
  end
  x, w
end

#if ! isinteractive()
#  if abspath(PROGRAM_FILE) == @__FILE__
#    setprecision(100) do         # standard double precision is 53
#      test_gauss_chebyshev2nd()
#    end
#  end
#end
