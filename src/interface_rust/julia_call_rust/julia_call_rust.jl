
const RUSTLIB = joinpath(@__DIR__, "target/debug/libmy_rust_lib")
println("Hello from Julia")
input = 10 #Int32(10)
output =  ccall(   #(:function or "function", "library"), Return type, (Input types,), arguments if any)
                (:double_input,
                RUSTLIB),
                Int32,          # Return type
                (Int32,),       # (Input types,)
                input)          # Arguments if any
println("As result of $input * 2 is: $output")
