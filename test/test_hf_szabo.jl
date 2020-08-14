
using SciAlgs: hartree_fock

@testset "Hartree-Fock: Szabo-Ostlund book" begin

   # TESTING H2
   total_energy, electronic_energy = hartree_fock([0., 1.4], [1, 1])
   szabo_energy = -1.8310
   @test isapprox(electronic_energy,szabo_energy,atol=1e-4)

   # TESTING HEH+
   total_energy, electronic_energy = hartree_fock([0., 1.4632], [2, 1])
   szabo_energy = -4.227529
   @test isapprox(electronic_energy,szabo_energy,atol=1e-6)
end

