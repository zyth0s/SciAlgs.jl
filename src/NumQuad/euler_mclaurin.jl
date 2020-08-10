
import SciAlgs.NumQuad: get_log10

# Calculates the nodes and weights of the Euler-McLaurin quadrature
# x_i = 
# w_i =
function euler_mclaurin(n)
  x = zeros(BigFloat,n)
  w = zeros(BigFloat,n)
  u = BigFloat(0)
  Δu =  BigFloat(1) / (n - 1)
  w = fill(Δu,n)
  w[1]   = Δu / 2
  w[end] = Δu / 2
  for i in 1:n
    x[i] = u
    u += Δu
  end
  x, w
end

#if ! isinteractive()
#  if abspath(PROGRAM_FILE) == @__FILE__
#    setprecision(100) do         # standard double precision is 53
#      test_euler_mclaurin()
#    end
#  end
#end
