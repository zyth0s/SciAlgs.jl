
# Must be out of if VERSION, otherwise => strange error
const FORTLIB    = joinpath(@__DIR__, "fortranlib.so")
const FORTSOURCE = joinpath(@__DIR__, "fortrancode.f90")

run(`gfortran -Wall -shared -O2 -fPIC -o $FORTLIB $FORTSOURCE`)

@static if VERSION >= v"1.5.0" # nicer syntax

   five = @ccall FORTLIB.__m_MOD_five()::Int
   println("Five = ",five)
else # works anyway

   println("Using function")
   five = ccall( (:__m_MOD_five, FORTLIB), Int, () )
   println("Five = ", five)
end
