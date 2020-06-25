
t = ccall(:clock, Int32, ())
println("Time: ",t)

four = ccall(:pow,Float64,(Float64,Float64),2,2)
println("Four = ",four)

