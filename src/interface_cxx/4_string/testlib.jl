# Load the module and generate the functions
module CppString
  using CxxWrap
  @wrapmodule(joinpath(@__DIR__, "build/lib/libtestlib"))

  function __init__()
    @initcxx
  end
end

@assert CppString.read_string("asd") == 3
