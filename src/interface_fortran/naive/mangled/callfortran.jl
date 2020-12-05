
# https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/#Passing-Pointers-for-Modifying-Inputs-1

# Must be out of if VERSION, otherwise => strange error
const FORTLIB    = joinpath(@__DIR__, "fortranlib.so")
const FORTSOURCE = joinpath(@__DIR__, "fortranlib.f90")

run(`gfortran -Wall -shared -O2 -fPIC -o $FORTLIB $FORTSOURCE`)

# Fortran passes all the arguments by reference, so we have to
# use references of the appropriate type.
v = Ref{Clong}(2)

# We are calling a subroutine, so return type is
# Void (for Fortran)  → Cvoid (via iso_c_binding) → Nothing (via Base.cconvert)
@static if VERSION >= v"1.5.0"  # clearer syntax

   @time @ccall FORTLIB.__basic_math_MOD_mypow(v::Ref{Clong})::Cvoid
else # works anyway

   # Argument types are specified with a tuple,
   # (Any,) isa Tuple => true
   # (Any)  isa Tuple => false
   ccall( (:__basic_math_MOD_mypow, FORTLIB), Cvoid, (Ref{Clong},), v )
end

# The content reference is accessed with []
@assert v[] == 4
println("Four = ",v[])
