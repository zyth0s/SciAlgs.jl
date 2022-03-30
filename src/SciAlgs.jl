module SciAlgs

include("car2int.jl")
include("geometry.jl")
include("io.jl")
include("validate_dni.jl")
include("validate_isbn.jl")

include("./Xtal/Xtal.jl")
include("./NumQuad/NumQuad.jl")
include("./Algorithms_from_THE_BOOK/TheBook.jl")
include("./ElStruct/McMurchieDavidson.jl")
include("./ElStruct/hf_szabo.jl")
include("./ElStruct/lattice_electrostatic_sum.jl")

end # module
