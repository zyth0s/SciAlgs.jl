
# References:
# [1] Clifford Boundary Conditions: A Simple Direct-Sum Evaluation of Madelung Constants
#     Nicolas Tavernier, Gian Luigi Bendazzoli, Vé ronique Brumas, Stefano Evangelisti,* and J. A. Berger*
# [2] Clifford boundary conditions for periodic systems: the Madelung constant of cubic crystals in 1, 2 and 3 dimensions
#     Nicolas Tavernier, Gian Luigi Bendazzoli, Véronique Brumas, Stefano Evangelisti & J. Arjan Berger*
# * ENVIRON (testing)

using Test
using Printf
using LinearAlgebra


@doc """
    madelung_clifford(abc, αβγ, positions, Q)

Compute the electrostatic sum of a distribution of point charges with Tavernier
et al. method (doi:10.1021/acs.jpclett.0c01684). The unit cell is defined by a
matrix containing the lattice vectors as rows. Ion coordinates are given as
columns of `positions`, and their charges are listed in `Q`.
The cell is assumed to have orthogonal axes (crystallographic cell).
"""
function madelung_clifford(L, positions, Q, Rmax=[120,120,120])

   # Lattice vectors as rows is the convention!!
   @assert size(L) == (3,3)       "Lattice-vectors matrix has wrong size"
   @assert size(positions,1) == 3 "Ion coordinates have to be stacked as columns"
   @assert size(positions,2) == length(Q) "Few/many charges or positions given"
   @assert isdiag(L'L) "Unit cell should have orthogonal axes"
   #@assert L[1,:] ⋅ L[2,:] ≈ L[2,:] ⋅ L[3,:] ≈ L[1,:] ⋅ L[3,:] ≈ 0 "Cell has no orthogonal axes"

   #G = L*L' # metric tensor
   #Ω = sqrt(det(G))

   Nat = length(Q)
   # Madelung constants
   M = zeros(Nat)

   # Sum over Clifford Super Cell (CSC), Equation (11)[1] or Equation (25)[2]
   for A in 1:Nat
      for B in 1:Nat

         # k1,k2,k3 in the paper
         for R1 in 0:Rmax[1]-1, # ∑'_R (': avoid self-interaction)
             R2 in 0:Rmax[2]-1,
             R3 in 0:Rmax[3]-1

            R = Float64[R1, R2, R3]

            R == [0,0,0] && A == B && continue

            # Squared distance between A and B
            Rᴬᴮ² = 0.0
            for d in 1:3

               Δx = positions[d,B] - positions[d,A] + R[d]

               # Squared equation (6)[1]
               Rᴬᴮ² += Rmax[d]^2 / (2π^2) * (1 - cos(2π/Rmax[d]*(Δx))) * norm(L[d,:])^2
               # Equation (2)[2]
               #Rᴬᴮ² += (Rmax[d] / π * sin(π*Δx/Rmax[d]))^2
            end
            M[A] += Q[B]/sqrt(Rᴬᴮ²)
         end
      end
   end

   # Equation (12)
   E = 0.5sum(Q .* M)
   #E /= cbrt(Ω) # if the sin formula for RAB is used
   @printf("E = %.15f\n",E)

   E
end

#   # Find the nearest-neighbor distance
#   R₀ = 99999.0
#   for A in 1:Nat
#      for B in 1:Nat
#
#         for R1 in -1:1, # ∑'_R (': avoid self-interaction)
#             R2 in -1:1,
#             R3 in -1:1
#
#            R = [R1, R2, R3]
#
#            R == [0,0,0] && A == B && continue
#
#            rij = positions[:,B] + R - positions[:,A]
#            R₀ = min(R₀, sqrt(rij'G*rij) )
#            #R₀ = min(R₀, norm(rij) )
#         end
#      end
#   end
#   @show R₀

# Test against https://github.com/lukeolson/pyewald/blob/master/pyewald/tests/test_ewaldsum.py

angs2bohr(angs) = inv(5.291_772_109_03e-1)*angs

@testset "Madelung Clifford: NaCl Madelung constant" begin

   # such that distance Na-Cl = 1 bohr.
   L = diagm([2,2,2]) # Rows are vectors
   positions = [0   0.5 0.5;
                0.5 0   0.5;
                0.5 0.5 0;
                0   0   0;
                0.5 0.5 0.5;
                0.5 0   0;
                0   0.5 0;
                0   0   0.5]'
   Q = [1, 1, 1, 1, -1, -1, -1, -1]

   E = madelung_clifford(L,positions,Q,[40,40,40])
   @test 0.25E ≈ −1.747_983_013_4 atol=1e-7 # from Table I K=40

   E = madelung_clifford(L,positions,Q,[60,60,60])
   @test 0.25E ≈ −1.747_750_568_2 atol=1e-7 # from Table I K=60

   E = madelung_clifford(L,positions,Q,[80,80,80])
   @test 0.25E ≈ −1.747_669_206_7 atol=1e-7 # from Table I K=80

   #@info "NaCl: Eelec = $E [Ha]"
   #@test 0.25E ≈ −1.747_564_594_633_182_190_636_212_035 atol=1e-4
   # Reference energy from Malik Mamode@DOI:10.1007/s10910-016-0705-9
end

@testset "Madelung Clifford: CsCl Madelung constant" begin

   L = diagm([1,1,1])
   positions = [0   0   0;
                0.5 0.5 0.5]'
   Q = [1, -1]

   E = madelung_clifford(L,positions,Q,[40,40,40])
   @test 0.5E * sqrt(3) ≈ −1.761_312_912_9 # Table 2 K=40

   E = madelung_clifford(L,positions,Q,[60,60,60])
   @test 0.5E * sqrt(3) ≈ −1.762_070_328_1 # Table 2 K=60

   #@debug "CsCl: Eelec = $E [Ha]"
   #@test 0.5E * sqrt(3) ≈ −1.762_674_773_070_99
end

@testset "Madelung Clifford: ZnS Madelung constant" begin

   acell = √2/2
   L = diagm([acell,acell,acell])
   positions = [0    0    0;
                0.25 0.25 0.25]'
   Q = [1, -1]

   # FIXME
   E = madelung_clifford(L,positions,Q)
   @test_broken 0.25E * sqrt(3) ≈ −1.638_066_314_9 # Table 3 K=40

   # FIXME
   E = madelung_clifford(L,positions,Q)
   @test_broken 0.25E * sqrt(3) ≈ −1.638_060_088_4 # Table 3 K=60

   #@debug "ZnS: Eelec = $E [Ha]"
   #@test 0.25E * sqrt(3) ≈ −1.638_055_053_388_79
end

@testset "Madelung Clifford: NaCl energy vs ENVIRON" begin

   acell = 10.944
   L = diagm([acell,acell,acell])
   positions = [0   0.5 0.5;    # Na
                0.5 0   0.5;    # Na
                0.5 0.5 0;      # Na
                0   0   0;      # Na
                0.5 0.5 0.5;    # Cl
                0.5 0   0;      # Cl
                0   0.5 0;      # Cl
                0   0   0.5]'   # Cl
   # QTAIM charges
   Q = [1, 1, 1, 1, -1, -1, -1, -1] .* 0.88

   E = madelung_clifford(L,positions,Q,[80,80,80])

   @debug "NaCl: Eelec = $E [Ha]"
   @test E ≈ -0.989264697 atol=1e-4 # Calculated with ENVIRON
end

@testset "Madelung Clifford: CsCl energy vs ENVIRON" begin

   acell = 7.7199
   L = diagm([acell,acell,acell])
   positions = [0   0   0;
                0.5 0.5 0.5]'
   # QTAIM charges
   Q = [1, -1] .* 0.83

   E = madelung_clifford(L,positions,Q,[80,80,80])

   @debug "CsCl: Eelec = $E [Ha]"
   @test E ≈ −0.181629364 atol=1e-4 # Calculated with ENVIRON
end

@testset "Madelung Clifford: LiCl energy vs ENVIRON" begin

   acell = 9.6933
   L = diagm([acell,acell,acell])
   positions = [0   0.5 0.5;    # Li
                0.5 0   0.5;    # Li
                0.5 0.5 0;      # Li
                0   0   0;      # Li
                0.5 0.5 0.5;    # Cl
                0.5 0   0;      # Cl
                0   0.5 0;      # Cl
                0   0   0.5]'   # Cl
   # QTAIM charges
   Q = [1, 1, 1, 1, -1, -1, -1, -1] .* 0.90

   E = madelung_clifford(L,positions,Q,[80,80,80])

   @debug "LiCl: Eelec = $E [Ha]"
   @test E ≈ -1.16825222 atol=1e-4 # Calculated with ENVIRON
end

@testset "Madelung Clifford: BN (cubic) energy vs ENVIRON" begin

   acell = 6.822
   L = diagm([acell,acell,acell])

   positions = [0   0.5 0.5;  # B; fractional coordinates
                0.5 0   0.5;  # B
                0.5 0.5 0;    # B
                0   0   0;    # B
                0.25 0.25 0.25;  # N
                0.25 0.75 0.75;  # N
                0.75 0.25 0.75;  # N
                0.75 0.75 0.25]' # N
   # QTAIM charges
   Q = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0] .* 2.16

   E = madelung_clifford(L,positions,Q,[80,80,80])

	@debug "BN (cubic): Eelec = $E [Ha]"
   @test E ≈ -10.3486490 atol=1e-4 # Calculated with ENVIRON
end

#@testset "Ewald: BN (hex.) energy vs ENVIRON" begin
#
#   acell = 4.7325
#   ccell = 12.5834
#   abc = [acell, acell, ccell]
#   αβγ = [90, 90, 120] .|> deg2rad
#
#   positions = [0.33333 0.66667 0.25;  # B; fractional coordinates
#                0.66667 0.33333 0.75;  # B
#                0.33333 0.66667 0.75;  # N
#                0.66667 0.33333 0.25]' # N
#   # QTAIM charges
#   Q = [1.0, 1.0, -1.0, -1.0] .* 2.214
#   ϵ = 1e-15
#
#   E = ewald(abc,αβγ,positions,Q,ϵ)
#
#	@debug "BN (hex.): Eelec = $E [Ha]"
#   @test E ≈ -5.53593310 atol=1e-6 # Calculated with ENVIRON
#end
#
#@testset "Ewald: MgB₂ energy vs ENVIRON" begin
#
#   acell = 5.8262
#   ccell = 6.6395
#   abc = [acell, acell, ccell]
#   αβγ = [90.0, 90.0, 120.0] .|> deg2rad
#   positions = [0       0       0;        # Mg
#                0.33333 0.66667 0.50000;  # B
#                0.66667 0.33333 0.50000]' # B
#   # QTAIM charges
#   Q = [2, -1, -1] * 0.81
#   ϵ = 1e-15
#
#   E = ewald(abc,αβγ,positions,Q,ϵ)
#
#	@debug "MgB₂: Eelec = $E [Ha]"
#   @test E ≈ -0.623919745 atol=1e-6 # Calculated with ENVIRON
#end
#
#
## Tested against values from the original paper
## doi:10.1103/PhysRevMaterials.2.013806
#
#@testset "Ewald: Al fcc vs CASTEP" begin
#
#   # Al; ICSD:43423
#   abc = [2.8636, 2.8636, 2.8636] .|> angs2bohr
#   αβγ = [60.0, 60.0, 60.0] .|> deg2rad
#
#   positions = zeros(3)
#   # Valence charge
#   Q = [3.0]
#   ϵ = 1e-15
#
#   E = ewald(abc,αβγ,positions,Q,ϵ)
#   @test E ≈ -2.695954572 # from Table I
#end
#
#@testset "Ewald: Si vs CASTEP" begin
#
#   # Si; ICSD:51688
#   abc = [3.8400, 3.8400, 3.8400] .|> angs2bohr
#   αβγ = [60.0, 60.0, 60.0] .|> deg2rad
#
#   positions = [0.00 0.00 0.00;  # Si
#                0.25 0.25 0.25]' # Si
#   # Valence charge
#   Q = [4.0, 4.0]
#   ϵ = 1e-15
#
#   E = ewald(abc,αβγ,positions,Q,ϵ)
#
#   @debug "Si: Eelec = $E [Ha]"
#   @test E ≈ -8.398574646 # from Table I
#end
#
#@testset "Ewald: SiO₂ vs CASTEP" begin
#
#   # SiO₂; ICSD:29122
#   abc = [4.9130, 4.9130, 5.4050] .|> angs2bohr
#   αβγ = [90.0, 90.0, 120.0] .|> deg2rad
#
#   positions = [0.41500  0.27200  0.21300; # O ; fractional coords
#                0.72800  0.14300  0.54633; # O
#                0.85700  0.58500  0.87967; # O
#                0.27200  0.41500  0.78700; # O
#                0.14300  0.72800  0.45367; # O
#                0.58500  0.85700  0.12033; # O
#                0.46500  0.00000  0.33333; # Si
#                0.00000  0.46500  0.66667; # Si
#                0.53500  0.53500  0.00000]'# Si
#   # Valence charges
#   Q = vcat(repeat([6.0], 6), repeat([4.0], 3))
#   ϵ = 1e-15
#
#   E = ewald(abc,αβγ,positions,Q,ϵ)
#
#   @debug "SiO₂: Eelec = $E [Ha]"
#   @test E ≈ -69.488098659 # from Table I
#end
#
#@testset "Electrostatic sum: Al₂SiO₅" begin
#
#   # Al₂SiO₅; ICSD:24275
#   L = diagm([14.7289033699982, 14.9260018049230, 10.5049875335275])
#
#   positions = [0.23030  0.13430  0.23900; # ; fractional coords
#                0.76970  0.86570  0.23900; #
#                0.26970  0.63430  0.26100; #
#                0.73030  0.36570  0.26100; #
#                0.76970  0.86570  0.76100; #
#                0.23030  0.13430  0.76100; #
#                0.73030  0.36570  0.73900; #
#                0.26970  0.63430  0.73900; #
#                0.00000  0.00000  0.24220; #
#                0.50000  0.50000  0.25780; #
#                0.00000  0.00000  0.75780; #
#                0.50000  0.50000  0.74220; #
#                0.37080  0.13870  0.50000; #
#                0.42320  0.36270  0.50000; #
#                0.62920  0.86130  0.50000; #
#                0.57680  0.63730  0.50000; #
#                0.12920  0.63870  0.00000; #
#                0.07680  0.86270  0.00000; #
#                0.87080  0.36130  0.00000; #
#                0.92320  0.13730  0.00000; #
#                0.24620  0.25290  0.00000; #
#                0.42400  0.36290  0.00000; #
#                0.10380  0.40130  0.00000; #
#                0.75380  0.74710  0.00000; #
#                0.57600  0.63710  0.00000; #
#                0.89620  0.59870  0.00000; #
#                0.25380  0.75290  0.50000; #
#                0.07600  0.86290  0.50000; #
#                0.39620  0.90130  0.50000; #
#                0.74620  0.24710  0.50000; #
#                0.92400  0.13710  0.50000; #
#                0.60380  0.09870  0.50000]'#
#   # Valence charges
#   Q = 6.0ones(size(positions, 2))    # O atoms (initial)
#   Q[[9,10,11,12,13,15,17,19]] .= 3.0 # Al atoms
#   Q[[21,24,27,30]]            .= 4.0 # Si atoms
#
#   E = madelung_clifford(L,positions,Q,[80,80,80])
#
#   @debug "Al₂SiO₅: Eelec = $E [Ha]"
#   @test E ≈ -244.055008450 atol=1e-4 # from Table I
#end

