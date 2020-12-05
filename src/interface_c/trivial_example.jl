
@static if VERSION >= v"1.5.0"

   # @ccall macro
   libc_abs(x::Int64  ) = @ccall abs(x::Cint   )::Cint
   libc_erf(x::Float64) = @ccall erf(x::Cdouble)::Cdouble

   @assert libc_abs(-4) == 4
else

   # ccall() function
   t = ccall(:clock, Int32, ())
   println("Time: ",t)

   four = ccall(:pow,Float64,(Float64,Float64),2,2)
   @assert four == 4
end
