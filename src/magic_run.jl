
# Run Julia files from the REPL with the following syntax:
# > @run "../file.jl" <-- String (any rel/abs path) [si autocompletion]
# > @run     file.jl  <-- Expr   (no ./ or ../ etc) [no autocompletion]
# > @run     file     <-- Symbol (no ./ or ../ etc) [no autocompletion]
macro run(filename)
   if typeof(filename) == String
      include(filename)
   elseif typeof(filename) == Expr
      include(join([string(filename.args[1]),
                    string(filename.args[2].value)],
                   "."))
   elseif typeof(filename) == Symbol
      include(string(filename) * ".jl")
   else
      error("The filepath must be a string.")
   end
end
