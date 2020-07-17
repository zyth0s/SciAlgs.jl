
# Taken from https://nbviewer.jupyter.org/github/mfherbst/course_julia_day/blob/master/07_Other_Language_Features.ipynb

code = """
double sum_array(double* array, int n) {
    double accu = 0.0;
    for (int i = 0; i < n; ++i) {
        accu += array[i];
    }
    return accu;
}
""";

open("sums.c", "w") do f
    write(f, code)
end
run(`cc -shared -o libsums.so sums.c`)
rm("sums.c")

v = [1.0, 2.0, 3.0]
res = ccall((:sum_array, "libsums.so"), Cdouble,
            (Ptr{Cdouble}, Cint), v, length(v))
