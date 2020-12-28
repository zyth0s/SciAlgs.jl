
# https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/#Passing-Pointers-for-Modifying-Inputs-1

# Must be out of if VERSION, otherwise => strange error
const FORTLIB    = joinpath(@__DIR__, "fortranlib.so")
const FORTSOURCE = joinpath(@__DIR__, "fortran_iso_c.f90")

run(`gfortran -Wall -shared -O2 -fPIC -o $FORTLIB $FORTSOURCE`)

# Fortran passes all the arguments by reference, so we have to
# use references of the appropriate type.
v = Ref{Clonglong}(2)

# We are calling a subroutine, so return type is
# Void (for Fortran)  → Cvoid (via iso_c_binding) → Nothing (via Base.cconvert)
@static if VERSION >= v"1.5.0"  # clearer syntax

   @ccall FORTLIB.mypow(v::Ref{Clonglong})::Cvoid
else # works anyway

   # Argument types are specified with a tuple,
   # (Any,) isa Tuple => true
   # (Any)  isa Tuple => false
   @time ccall( (:mypow, FORTLIB), Cvoid, (Ref{Clonglong},), v )
end

# The content reference is accessed with []
@assert v[] == 4
println("Four = ",v[])


# TODO give examples of all combinations possible:
#      * subroutine call passing integer by reference (inout)
#      * subroutine call passing array   by reference (inout)
#      * subroutine call passing integer by reference (in)
#      * subroutine call passing array   by reference (in)
#      * function call that returns integer
#      * function call that returns array
