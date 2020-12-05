
# Taken from https://nbviewer.jupyter.org/github/mfherbst/course_julia_day/blob/master/07_Other_Language_Features.ipynb

# C code string
code = """
double sum_array(double* array, int n) {
    double accu = 0.0;
    for (int i = 0; i < n; ++i) {
        accu += array[i];
    }
    return accu;
}
""";

# save it to a file
open("sums.c", "w") do f
    write(f, code)
end
# compile it
run(`cc -shared -o libsums.so sums.c`)
rm("sums.c")

# Test

const libsums = "libsums.so"

@static if VERSION >= v"1.5.0"

   v = [1.0, 2.0, 3.0]

   res = @ccall libsums.sum_array(v::Ptr{Cdouble}, length(v)::Cint)::Cdouble

   @assert res ≈ 6
else

   v = [1.0, 2.0, 3.0]

   res = ccall((:sum_array, libsums), Cdouble,
               (Ptr{Cdouble}, Cint), v, length(v))
   @assert res ≈ 6
end
