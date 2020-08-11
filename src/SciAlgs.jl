module SciAlgs

include("magic_run.jl")

include("./Xtal/Xtal.jl")
include("./NumQuad/NumQuad.jl")
include("./Algorithms_from_THE_BOOK/TheBook.jl")

greet() = print("Hello World!daniel")

#SciAlgs.Xtal.example1()

end # module
