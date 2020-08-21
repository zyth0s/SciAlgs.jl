
#using Libdl
#push!(Libdl.DL_LOAD_PATH,pwd())

@static if VERSION >= v"1.5.0"

   const libmean = "libmean.so"
   g(x::Float64, y::Float64) = @ccall libmean.mean(x::Cdouble, y::Cdouble)::Cdouble

   @assert g(4.0,6.0) ≈ 5
else

   f(x::Float64,y::Float64) = ccall((:mean,"libmean.so"),Float64,(Float64,Float64),x,y)
   @assert f(4.0,6.0) ≈ 5
end
