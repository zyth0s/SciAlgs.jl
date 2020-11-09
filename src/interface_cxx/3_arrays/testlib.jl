# Load the module and generate the functions
module CppArrays
  using CxxWrap
  @wrapmodule(joinpath(@__DIR__, "build/lib/libtestlib"))

  function __init__()
    @initcxx
  end
end
using .CppArrays

# Take Julia array as argument
ta = [1.0, 2.0]
CppArrays.test_array_set(ta, 0, 3.0)
@assert ta ≈ [3.0, 2.0]

# Produce a 1D Julia array

# Produce a nD Julia array
#display(CppClass.const_matrix())
