
#using Libdl
#push!(Libdl.DL_LOAD_PATH,pwd())

# Must be out of if VERSION, otherwise => strange error
const CLIB    = joinpath(@__DIR__, "libmean.so")
const CSOURCE = joinpath(@__DIR__, "calc_mean.c")

run(`gcc -Wall -shared -O2 -fPIC -o $CLIB $CSOURCE`)

@static if VERSION >= v"1.5.0"

   g(x::Float64, y::Float64) = @ccall CLIB.mean(x::Cdouble, y::Cdouble)::Cdouble

   @time g(4.0,6.0)
   @assert g(4.0,6.0) ≈ 5
else

   f(x::Float64,y::Float64) = ccall((:mean, CLIB),Float64,(Float64,Float64),x,y)
   @assert f(4.0,6.0) ≈ 5
end
