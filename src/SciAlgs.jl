module SciAlgs

# Run Julia files from the REPL e.g. > @run file.jl
macro run(filename)
   if typeof(filename) == Expr
      include(repr(filename.args[1])[2:end] * ".jl")
   elseif typeof(filename) == String
      include(filename)
   end
end

module Bar
   greet2() = print("Hello World!")
end
include("./Xtal/Xtal.jl")

greet() = print("Hello World!daniel")

#SciAlgs.Xtal.example1()

include("extra_file.jl")

end # module
