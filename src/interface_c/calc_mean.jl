
#using Libdl
#push!(Libdl.DL_LOAD_PATH,pwd())

f(x::Float64,y::Float64) = ccall((:mean,"libmean.so"),Float64,(Float64,Float64),x,y)

@show f(4.0,6.0)
