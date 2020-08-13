
# References:
# - (HJO) "Molecular Electronic Structure Theory" 1st ed. by T. Helgaker, P. Jorgensen, J. Olsen
#         is the default reference
# - (FV) "Fundamentals of Molecular Integrals Evaluation" by J. T. Fermann, E. Valeev
# Note: This implementation is influenced by Joshua Goings Python implementation

module BasisSet

export BasisFunction, normalize_basis!


import Combinatorics: doublefactorial

abstract type BasisFunctionType end

mutable struct BasisFunction <: BasisFunctionType
   origin; shell; exps; coefs; norm
end

function normalize_basis!(basis::BasisFunction)
   # HJO gives no single expression, see (9.2.4), (8.2.14) 
   l,m,n = basis.shell
   L = l + m + n
   # Primitive (PGBFs) normalization ∼ (2.11) in FV
   basis.norm = sqrt.( 2^(2*L+1.5)*basis.exps.^(L+1.5)/
                      dfact2(2l-1)/dfact2(2m-1)/dfact2(2n-1)/
                      (π^1.5))

   # Contracted (CGBFs) normalization (for s functions first)
   # (2.25), also (2.18) of FV
   prefactor = π^1.5 *dfact2(2l-1)*dfact2(2m-1)*dfact2(2n-1)/2^L
   N = 0.0
   for ia in eachindex(basis.exps), ib in eachindex(basis.exps)
       # for contracted functions same ang. momentum!! (2.25) in FV
       N += basis.norm[ia]*basis.norm[ib]*basis.coefs[ia]*basis.coefs[ib]/
       (basis.exps[ia]+basis.exps[ib])^(L+1.5)
   end
   N *= prefactor
   N = 1/√N
   basis.coefs *= N
   nothing
end

dfact2(x::Int) = (x == -1) ? 1.0 : Float64(doublefactorial(x)) # BigInt -> Float64


end # module
