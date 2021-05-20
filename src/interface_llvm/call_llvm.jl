
# https://llvm.org/docs/LangRef.html#llvm-exp-intrinsic
llvm_exp(x) = ccall("llvm.exp.f64", llvmcall, Float64, (Float64,), x);
llvm_exp(3.0)

# https://llvm.org/docs/LangRef.html#llvm-pow-intrinsic
llvm_pow(a, b) = ccall("llvm.pow.f64", llvmcall, Float64, (Float64, Float64), a, b);
llvm_pow(3.0, 4.0)

# http://kristofferc.github.io/post/intrinsics/
const __m128i = NTuple{2, VecElement{Int64}};
aesdec(a, roundkey) = ccall("llvm.x86.aesni.aesdec", llvmcall, __m128i, (__m128i, __m128i), a, roundkey);
aesdec(__m128i((213132, 13131)), __m128i((31231, 43213)))

