
# Run Julia files from the REPL, e.g.
# julia> @run file.jl      (only files in same folder)
# julia> @run "../file.jl" (has autocompletion)
macro run(filename)
   if typeof(filename) == Expr
      include(repr(filename.args[1])[2:end] * ".jl")
   elseif typeof(filename) == String
      include(filename)
   end
end
