# Load the module and generate the functions
module CppClass
  using CxxWrap
  @wrapmodule(joinpath(@__DIR__, "build/lib/libtestlib"))

  function __init__()
    @initcxx
  end
end

# Test
f = CppClass.Foo(4)
@assert CppClass.add(f, 4) == 8 # 4 + 4
