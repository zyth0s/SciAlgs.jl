# Load the module and generate the functions
module CppHello
  using CxxWrap
  @wrapmodule(joinpath(@__DIR__, "build/lib/libtestlib"))

  function __init__()
    @initcxx
  end
end

# Call greet and show the result
@assert CppHello.greet() |> unsafe_string == "Hello"
