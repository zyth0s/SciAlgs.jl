
# References:
# * https://www.nist.gov/mml/csd/chemical-informatics-group/spce-water-reference-calculations-non-cuboid-cell-10a-cutoff
#       TODO H₂O tests
# * 5.5.2 in M. P. Allen and D. J. Tildesley, Computer Simulation of Liquids (Oxford University Press, New York, 1989).
# * https://ccse.jaea.go.jp/software/PIMD/doc/manual/node575.html
# * ENVIRON (testing)
# * Cutoffs as in qewald.m Copyright (C) 2009-2011 A. Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
#                                              and V. Lua~na <victor@carbono.quimica.uniovi.es>. Universidad de Oviedo.

using Test
using Printf
using LinearAlgebra

erfc(x::Float64) = @ccall erfc(x::Float64)::Float64

function cutoffs_r(abc, αβγ, Nat, ∑Q², α, Ω, ϵ; SGROW = 1.4, EPSCUT = 1e-5)

   r_cut1 = 1
   r_cut2 = 2/ SGROW

   while true
      r_cut2 *= SGROW
      error = π * Nat^2 * ∑Q² * erfc(α * r_cut2) / Ω / α^2
      error < ϵ && break
   end

   while (r_cut2-r_cut1) > EPSCUT
      r_cut = 0.5(r_cut1+r_cut2)
      error = π * Nat^2 * ∑Q² * erfc(α * r_cut) / Ω / α^2
      if error > ϵ
        r_cut1 = r_cut
      else
        r_cut2 = r_cut
      end
   end
   r_cut = 0.5(r_cut1+r_cut2);
   #@printf("r_cut = %.15f\n",r_cut)

   # real space cells to explore
   Rmax = zeros(3)
   a, b, c = abc;  sinα, sinβ, sinγ = sin.(αβγ)
   Rmax[1] = b * c * sinα
   Rmax[2] = a * c * sinβ
   Rmax[3] = a * b * sinγ
   Rmax = floor.(r_cut * Rmax / Ω) .+ 1
   #@printf("Rmax = %d %d %d\n",Rmax...)

   r_cut, Rmax
end

function cutoffs_h(abc, Nat, ∑Q², α, ϵ; SGROW = 1.4, EPSCUT = 1e-5)

   h_cut1 = 1
   h_cut2 = 2 / SGROW
   while true
      h_cut2 *= SGROW
      error = Nat^2 * ∑Q² * α * erfc(0.5h_cut2 / α) / √π
      error < ϵ && break
   end
   while (h_cut2-h_cut1) > EPSCUT
      h_cut = 0.5(h_cut1+h_cut2)
      error = Nat^2 * ∑Q² * α * erfc(0.5h_cut / α) / √π
      if error > ϵ
         h_cut1 = h_cut
      else
         h_cut2 = h_cut
      end
   end
   h_cut = 0.5(h_cut1+h_cut2)
   #@printf("h_cut = %.15f\n",h_cut)

   # reciprocal space cells to explore
   Kmax = floor.(abc .* h_cut ./ (2π)) .+ 1
   #@printf("Kmax = %d %d %d\n",Kmax...)

   h_cut, Kmax
end

@doc """
    ewald(abc, αβγ, positions, Q, ϵ, α=missing)

Compute the electrostatic sum of a distribution of point charges with
Ewald's method. The unit cell is defined by the cell parameters, lengths `abc`
and angles `αβγ`. Ion coordinates are given as columns of `positions`, and
their charges are listed in `Q`. The parameter `ϵ` controls the admitted error
and `α` is the short/long range split parameter.
"""
function ewald(abc, αβγ, positions, Q, ϵ=10^(log10(eps(Float64))+4), α=missing)

   @assert size(abc) == (3,)      "abc has to be a column vector"
   @assert size(αβγ) == (3,)      "αβγ has to be a column vector"
   @assert size(positions,1) == 3 "Ion coordinates have to be stacked as columns"
   @assert size(positions,2) == length(Q) "Few/many charges or positions given"

   G = begin
      a, b, c = abc;  cosα, cosβ, cosγ = cos.(αβγ)
      [ a*a       a*b*cosγ  a*c*cosβ;
        a*b*cosγ  b*b       b*c*cosα;
        a*c*cosβ  b*c*cosα  c*c     ]
   end
   #@info "Metric tensor G" G

   Ω = √(det(G))
   Ginv = inv(G)
   Nat = length(Q)
   ∑Q  = sum(Q)
   ∑Q² = sum(Q.^2)
   #!(isapprox(∑Q, 0, atol=1e-8)) && @warn "Requires total charge neutrality; ∑Q = $∑Q"

   # calculate split parameter, short/long range
   α = ismissing(α) ? sqrt(π*abc[2]*sin(αβγ[3])/Ω) : α
   #@printf("α = %.15f\n",α)

   # real & reciprocal space cutoffs; based on the decay of exp() & erfc()
   r_cut, Rmax = cutoffs_r(abc, αβγ, Nat, ∑Q², α, Ω, ϵ)
   h_cut, Kmax = cutoffs_h(abc,      Nat, ∑Q², α,    ϵ)

   # Real-space sum
   ∑Eᵣ = 0
   for i in 1:Nat,
       j in 1:Nat

      for R1 in -Rmax[1]:Rmax[1], # ∑'_R (': avoid self-interaction)
          R2 in -Rmax[2]:Rmax[2],
          R3 in -Rmax[3]:Rmax[3]

         R = [R1, R2, R3]

         # i is in the central cell (R=(0,0,0)) whereas j is at any cell (rⱼ + R)
         𝐫ᵢⱼ = positions[:, i] - positions[:,j] - R
         # Distance |𝐫ᵢⱼ|²
         rᵢⱼ = sqrt(𝐫ᵢⱼ'G*𝐫ᵢⱼ)
         # Discard self-interaction (R=[0,0,0] && i==j) and far ions
         1e-12 < rᵢⱼ < r_cut || continue
         ∑Eᵣ += Q[i] * Q[j] * erfc(α*rᵢⱼ) / rᵢⱼ
      end
   end
   ∑Eᵣ /= 2 # correct double counting
   #@printf("∑Eᵣ = %.15f\n",∑Eᵣ)

   # Reciprocal-space sum
   ∑Eₖ = 0
   for K1 in -Kmax[1]:Kmax[1], # ∑_(𝐊≠(0,0,0))
       K2 in -Kmax[2]:Kmax[2],
       K3 in -Kmax[3]:Kmax[3]

      h = 2π * [K1, K2, K3] # = 2π 𝐤 : crystallographic vector
      # Norm of reciprocal vector, |𝐡|
      h_norm = sqrt(h'Ginv*h)
      1e-12 < h_norm < h_cut || continue

      # Structure factor: 𝐒(𝐡) = ∑ᵢ qᵢ exp(-i 𝐡 ⋅ 𝐫ᵢ)
      Sh = 0
      for i in 1:Nat
         Sh += Q[i] * exp(-im*h'positions[:,i])
      end
      Sh² = Sh'Sh # S(h) ⋅ S(-h) = |S(h)|²

      exponent = (0.5h_norm / α)^2 # = (π k / α)^2

      ∑Eₖ += Sh² / h_norm^2 * exp(-exponent)
   end
   ∑Eₖ *= 2π / Ω
   #@printf("∑Eₖ = %.15f\n",∑Eₖ)

   # 𝐤 = 0, self-energy term
   ∑₀ = - α * ∑Q² / √π
   #@printf("∑₀ = %.15f\n",∑₀)

   # compensating background charge term for charge neutrality
   ∑E_back = - 0.5π * ∑Q^2 / α^2 / Ω

   ∑Eₖ + ∑Eᵣ + ∑₀ + ∑E_back
end

# Testing
# =======

angs2bohr(angs) = inv(5.291_772_109_03e-1)*angs

@testset "Ewald: NaCl Madelung constant" begin

   acell = 2 # such that distance Na-Cl = 1 bohr.
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad
   positions = [0   0.5 0.5;
                0.5 0   0.5;
                0.5 0.5 0;
                0   0   0;
                0.5 0.5 0.5;
                0.5 0   0;
                0   0.5 0;
                0   0   0.5]'
   Q = [1, 1, 1, 1, -1, -1, -1, -1]

   E = ewald(abc,αβγ,positions,Q)

   @debug "NaCl: Eelec = $E [Ha]"
   @test 0.25E ≈ −1.747_564_594_633_182_190_636_212_035
   # Reference energy from Malik Mamode@DOI:10.1007/s10910-016-0705-9
   # also oeis.org/A085469
end

@testset "Ewald: CsCl Madelung constant" begin

   acell = 1
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad
   positions = [0   0   0;
                0.5 0.5 0.5]'
   Q = [1, -1]

   E = ewald(abc,αβγ,positions,Q)

   @debug "CsCl: Eelec = $E [Ha]"
   @test 0.5E * sqrt(3) ≈ −1.762_674_773_070_99
   # http://oeis.org/A181152
end

@testset "Ewald: ZnS Madelung constant" begin

   acell = √2/2
   abc = [acell, acell, acell]
   αβγ = [60, 60, 60] .|> deg2rad
   positions = [0    0    0;
                0.25 0.25 0.25]'
   Q = [1, -1]

   E = ewald(abc,αβγ,positions,Q)

   @debug "ZnS: Eelec = $E [Ha]"
   @test 0.25E * sqrt(3) ≈ −1.638_055_053_388_79
   # http://oeis.org/A182566
end

@testset "Ewald: NaCl energy vs ENVIRON" begin

   acell = 10.944
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad
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

   E = ewald(abc,αβγ,positions,Q)

   @debug "NaCl: Eelec = $E [Ha]"
   @test E ≈ -0.989264697 atol=1e-6 # Calculated with ENVIRON
end

@testset "Ewald: CsCl energy vs ENVIRON" begin

   acell = 7.7199
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad
   positions = [0   0   0;
                0.5 0.5 0.5]'
   # QTAIM charges
   Q = [1, -1] .* 0.83

   E = ewald(abc,αβγ,positions,Q)

   @debug "CsCl: Eelec = $E [Ha]"
   @test E ≈ −0.181629364 atol=1e-6 # Calculated with ENVIRON
end

@testset "Ewald: LiCl energy vs ENVIRON" begin

   acell = 9.6933
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad
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

   E = ewald(abc,αβγ,positions,Q)

   @debug "LiCl: Eelec = $E [Ha]"
   @test E ≈ -1.16825222 atol=1e-6 # Calculated with ENVIRON
end

@testset "Ewald: BN (cubic) energy vs ENVIRON" begin

   acell = 6.822
   abc = [acell, acell, acell]
   αβγ = [90, 90, 90] .|> deg2rad

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

   E = ewald(abc,αβγ,positions,Q)

	@debug "BN (cubic): Eelec = $E [Ha]"
   @test E ≈ -10.3486490 atol=1e-6 # Calculated with ENVIRON
end

@testset "Ewald: BN (hex.) energy vs ENVIRON" begin

   acell = 4.7325
   ccell = 12.5834
   abc = [acell, acell, ccell]
   αβγ = [90, 90, 120] .|> deg2rad

   positions = [0.33333 0.66667 0.25;  # B; fractional coordinates
                0.66667 0.33333 0.75;  # B
                0.33333 0.66667 0.75;  # N
                0.66667 0.33333 0.25]' # N
   # QTAIM charges
   Q = [1.0, 1.0, -1.0, -1.0] .* 2.214

   E = ewald(abc,αβγ,positions,Q)

	@debug "BN (hex.): Eelec = $E [Ha]"
   @test E ≈ -5.53593310 atol=1e-6 # Calculated with ENVIRON
end

@testset "Ewald: MgB₂ energy vs ENVIRON" begin

   acell = 5.8262
   ccell = 6.6395
   abc = [acell, acell, ccell]
   αβγ = [90.0, 90.0, 120.0] .|> deg2rad
   positions = [0       0       0;        # Mg
                0.33333 0.66667 0.50000;  # B
                0.66667 0.33333 0.50000]' # B
   # QTAIM charges
   Q = [2, -1, -1] * 0.81

   E = ewald(abc,αβγ,positions,Q)

	@debug "MgB₂: Eelec = $E [Ha]"
   @test E ≈ -0.623919745 atol=1e-6 # Calculated with ENVIRON
end


# Tested against values from the original paper
# doi:10.1103/PhysRevMaterials.2.013806

@testset "Ewald: Al fcc vs CASTEP" begin

   # Al; ICSD:43423
   abc = [2.8636, 2.8636, 2.8636] .|> angs2bohr
   αβγ = [60.0, 60.0, 60.0] .|> deg2rad

   positions = zeros(3)
   # Valence charge
   Q = [3.0]

   E = ewald(abc,αβγ,positions,Q)
   @test E ≈ -2.695954572 # from Table I
end

@testset "Ewald: Si vs CASTEP" begin

   # Si; ICSD:51688
   abc = [3.8400, 3.8400, 3.8400] .|> angs2bohr
   αβγ = [60.0, 60.0, 60.0] .|> deg2rad

   positions = [0.00 0.00 0.00;  # Si
                0.25 0.25 0.25]' # Si
   # Valence charge
   Q = [4.0, 4.0]

   E = ewald(abc,αβγ,positions,Q)

   @debug "Si: Eelec = $E [Ha]"
   @test E ≈ -8.398574646 # from Table I
end

@testset "Ewald: SiO₂ vs CASTEP" begin

   # SiO₂; ICSD:29122
   abc = [4.9130, 4.9130, 5.4050] .|> angs2bohr
   αβγ = [90.0, 90.0, 120.0] .|> deg2rad

   positions = [0.41500  0.27200  0.21300; # O ; fractional coords
                0.72800  0.14300  0.54633; # O
                0.85700  0.58500  0.87967; # O
                0.27200  0.41500  0.78700; # O
                0.14300  0.72800  0.45367; # O
                0.58500  0.85700  0.12033; # O
                0.46500  0.00000  0.33333; # Si
                0.00000  0.46500  0.66667; # Si
                0.53500  0.53500  0.00000]'# Si
   # Valence charges
   Q = vcat(repeat([6.0], 6), repeat([4.0], 3))

   E = ewald(abc,αβγ,positions,Q)

   @debug "SiO₂: Eelec = $E [Ha]"
   @test E ≈ -69.488098659 # from Table I
end

@testset "Electrostatic sum: Al₂SiO₅" begin

   # Al₂SiO₅; ICSD:24275
   abc = [14.7289033699982, 14.9260018049230, 10.5049875335275]
   #abc = [7.7942, 7.8985, 5.5590] .|> angs2bohr
   αβγ = [90.0, 90.0, 90.0] .|> deg2rad

   positions = [0.23030  0.13430  0.23900; # ; fractional coords
                0.76970  0.86570  0.23900; #
                0.26970  0.63430  0.26100; #
                0.73030  0.36570  0.26100; #
                0.76970  0.86570  0.76100; #
                0.23030  0.13430  0.76100; #
                0.73030  0.36570  0.73900; #
                0.26970  0.63430  0.73900; #
                0.00000  0.00000  0.24220; #
                0.50000  0.50000  0.25780; #
                0.00000  0.00000  0.75780; #
                0.50000  0.50000  0.74220; #
                0.37080  0.13870  0.50000; #
                0.42320  0.36270  0.50000; #
                0.62920  0.86130  0.50000; #
                0.57680  0.63730  0.50000; #
                0.12920  0.63870  0.00000; #
                0.07680  0.86270  0.00000; #
                0.87080  0.36130  0.00000; #
                0.92320  0.13730  0.00000; #
                0.24620  0.25290  0.00000; #
                0.42400  0.36290  0.00000; #
                0.10380  0.40130  0.00000; #
                0.75380  0.74710  0.00000; #
                0.57600  0.63710  0.00000; #
                0.89620  0.59870  0.00000; #
                0.25380  0.75290  0.50000; #
                0.07600  0.86290  0.50000; #
                0.39620  0.90130  0.50000; #
                0.74620  0.24710  0.50000; #
                0.92400  0.13710  0.50000; #
                0.60380  0.09870  0.50000]'#
   # Valence charges
   Q = 6.0ones(size(positions, 2))    # O atoms (initial)
   Q[[9,10,11,12,13,15,17,19]] .= 3.0 # Al atoms
   Q[[21,24,27,30]]            .= 4.0 # Si atoms

   E = ewald(abc,αβγ,positions,Q)

   @debug "Al₂SiO₅: Eelec = $E [Ha]"
   @test E ≈ -244.055008450 # from Table I
end

