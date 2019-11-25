
module NumQuad
   
# Newton-Cotes family
include("trapezoidal.jl")
include("euler_mclaurin.jl")
# Clenshaw-Curtis family
include("clenshaw_curtis.jl")
# Gauss family
include("gauss_legendre.jl")
include("gauss_chebyshev1st.jl")
include("gauss_chebyshev2nd.jl")
include("test_laguerre.jl")
include("test_chebyshev_integration.jl")
# Non-classical Gauss
include("maxwell.jl")
include("multiexp.jl")
# Other
include("perez_jorda.jl")

end
