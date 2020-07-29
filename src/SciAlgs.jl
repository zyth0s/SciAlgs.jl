module SciAlgs

include("magic_run.jl")

include("./Xtal/Xtal.jl")
include("./NumQuad/NumQuad.jl")
using SciAlgs.NumQuad: gauss_legendre

greet() = print("Hello World!daniel")

#SciAlgs.Xtal.example1()

include("extra_file.jl")

end # module
