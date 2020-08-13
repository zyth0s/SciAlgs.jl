
using SciAlgs.McMurchieDavidson: BasisFunction, normalize_basis!,
                                 S, T, V, ERI

@testset "McMurchie-Davidson scheme: s-function" begin

   origin = [1.0, 2.0, 3.0]
   shell  = (0,0,0) # e.g. pₓ-orbitals are (1,0,0)
   exps   = [3.42525091, 0.62391373, 0.16885540]
   coefs  = [0.15432897, 0.53532814, 0.44463454]

   a = BasisFunction(origin,shell,exps,coefs,missing)
   normalize_basis!(a)

   @test S(a,a)               ≈ 1.0                #|| error("Wrong overlap")
   @test T(a,a)               ≈ 0.760031883566609  #|| error("Wrong kinetic energy")
   @test V(a,a,[1.0,0.0,0.0]) ≈ 0.2771825991512926 #|| error("Wrong rep energy")
   @test ERI(a,a,a,a)         ≈ 0.7746059439198977 #|| error("Wrong eri")
end

