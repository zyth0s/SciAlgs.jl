
const RUSTLIB = joinpath(@__DIR__, "target/release/libmy_rust_lib")
println("Hello from Julia!")
input = 10 #Int32(10)
output =  ccall(   #(:function or "function", "library"), Return type, (Input types,), arguments if any)
                (:double_input,
                RUSTLIB),
                Int32,          # Return type
                (Int32,),       # (Input types,)
                input)          # Arguments if any
println()
println("The result of $input * 2 is $output.")

println()
v = zeros(3)
@ccall RUSTLIB.modify_vectorf64(v::Ptr{Float64},
                                length(v)::Int64)::Cvoid
println("Vector taken from Rust: v = $v")


println()
m = zeros(3,5)
@ccall RUSTLIB.modify_matrixf64(m::Ptr{Float64},
                                size(m,1)::Int64,
                                size(m,2)::Int64)::Cvoid
@info "Matrix taken from Rust:" m
#@info "Matrix taken from Rust: m' = (transposed for comparison)" m'
