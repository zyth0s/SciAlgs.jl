
# Must be out of if VERSION, otherwise => strange error
const FORTLIB = "fortranlib.so"

@static if VERSION >= v"1.5.0" # nicer syntax

   five = @ccall FORTLIB.__m_MOD_five()::Int
   println("Five = ",five)
else # works anyway

   println("Using function")
   five = ccall( (:__m_MOD_five, FORTLIB), Int, () )
   println("Five = ", five)
end
