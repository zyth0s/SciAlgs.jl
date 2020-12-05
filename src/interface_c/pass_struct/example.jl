C_header = """
struct AIMWFN_Mol_Gamess_HF {
   int nmo;
};
int foo(struct AIMWFN_Mol_Gamess_HF wfn);
"""
C_code = """
#include <stdio.h>
#include "struct.h"
int foo(struct AIMWFN_Mol_Gamess_HF wfn) {
  //printf("%5d\\n", wfn.nmo);
  return wfn.nmo;
}
"""
open("struct.h", "w") do f
    write(f, C_header)
end
open("struct.c", "w") do f
    write(f, C_code)
end
const LIBSTRUCT = joinpath(@__DIR__, "libstruct.so")
run(`gcc -Wall -shared -O2 -fPIC -o $LIBSTRUCT struct.c`)
rm("struct.c")
rm("struct.h")

# ------------------------------------------------------------
# Automatically generated using Clang.jl
struct AIMWFN_Mol_Gamess_HF
    nmo::Cint
end
function foo(wfn)
    ccall((:foo, LIBSTRUCT), Cint, (AIMWFN_Mol_Gamess_HF,), wfn)
end
# ------------------------------------------------------------

nmo = 8
wfn = AIMWFN_Mol_Gamess_HF(nmo)
@assert foo(wfn) == nmo
