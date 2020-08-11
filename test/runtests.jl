using Test

# NumQuad
include("test_gauss_legendre.jl")
include("test_clenshaw_curtis.jl")
include("test_gauss_chebyshev1st.jl")
include("test_gauss_chebyshev2nd.jl")
include("test_euler_mclaurin.jl")
include("test_trapezoidal.jl")

# Algorithms from THE BOOK
include("test_archimedes.jl")
