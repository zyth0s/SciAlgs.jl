
using SciAlgs: lattice_electrostatic_sum

# Tested against values from the original paper
# doi:10.1103/PhysRevMaterials.2.013806

@testset "Electrostatic sum: Al fcc" begin

   Rd_hat = 2.0 # parameter from Table I

   # Al; ICSD:43423
   a₁ = [5.41141973394663 , 0.0              , 0.0]
   a₂ = [2.70570986697332 , 4.68642696013821 , 0.0]
   a₃ = [2.70570986697332 , 1.56214232004608 , 4.41840571073226]

   r_ions = zeros(1,3)
   # Valence charge
   Qval = [3.0]

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @test E ≈ -2.695954572 # from Table I
end

@testset "Electrostatic sum: Si" begin

   Rd_hat = 2.0 # parameter from Table I

   # Si; ICSD:51688
   a₁ = [7.25654832321381, 0.00000000000000, 0.00000000000000]
   a₂ = [3.62827416160690, 6.28435519169252, 0.00000000000000]
   a₃ = [3.62827416160690, 2.09478506389751, 5.92494689524090]

   r_ions = [0.00  0.00  0.00; # Si ; fractional coords
             0.25  0.25  0.25] # Si
   #r_ions = (hcat(a₁,a₂,a₃)*r_ions'ᵀ)'ᵀ # cartesian
   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Valence charges
   Qval = [4.0, 4.0]

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @debug "Si: Eelec = $E [Ha]"
   @test E ≈ -8.398574646 # from Table I
end

@testset "Electrostatic sum: SiO₂" begin

   Rd_hat = 2.0 # parameter from Table I

   # SiO₂; ICSD:29122
   a₁ = [ 9.28422445623683, 0.00000000000000, 0.00000000000000]
   a₂ = [-4.64211222811842, 8.04037423353787, 0.00000000000000]
   a₃ = [ 0.00000000000000, 0.00000000000000, 10.2139697101486]

   r_ions = [0.41500  0.27200  0.21300; # O ; fractional coords
             0.72800  0.14300  0.54633; # O
             0.85700  0.58500  0.87967; # O
             0.27200  0.41500  0.78700; # O
             0.14300  0.72800  0.45367; # O
             0.58500  0.85700  0.12033; # O
             0.46500  0.00000  0.33333; # Si
             0.00000  0.46500  0.66667; # Si
             0.53500  0.53500  0.00000] # Si
   #r_ions = (hcat(a₁,a₂,a₃)*r_ions'ᵀ)'ᵀ # cartesian
   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Valence charges
   Qval = vcat(repeat([6.0], 6), repeat([4.0], 3))

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @debug "SiO₂: Eelec = $E [Ha]"
   @test E ≈ -69.488098659 # from Table I
end

@testset "Electrostatic sum: Al₂SiO₅" begin

   Rd_hat = 2.0 # parameter from Table I

   # Al₂SiO₅; ICSD:24275
   a₁ = [14.7289033699982, 0.00000000000000, 0.00000000000000]
   a₂ = [0.00000000000000, 14.9260018049230, 0.00000000000000]
   a₃ = [0.00000000000000, 0.00000000000000, 10.5049875335275]

   r_ions = [0.23030  0.13430  0.23900; # ; fractional coords
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
             0.60380  0.09870  0.50000] #
   #r_ions = (hcat(a₁,a₂,a₃)*r_ions'ᵀ)'ᵀ # cartesian
   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Valence charges
   Qval = 6.0ones(size(r_ions, 1)) # O atoms (initial)
   Qval[[9,10,11,12,13,15,17,19]] .= 3.0 # Al atoms
   Qval[[21,24,27,30]]            .= 4.0 # Si atoms

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @debug "Al₂SiO₅: Eelec = $E [Ha]"
   @test E ≈ -244.055008450 # from Table I
end

@testset "Electrostatic sum: NaCl" begin

   Rd_hat = 2.0 # parameter from Table I

   # NaCl; d(Na-Cl) = 1.0 bohr
   a₁ = [2.0, 0.0, 0.0]
   a₂ = [0.0, 2.0, 0.0]
   a₃ = [0.0, 0.0, 2.0]

   r_ions = [0   0.5 0.5;  # Na; fractional coordinates
             0.5 0   0.5;  # Na
             0.5 0.5 0;    # Na
             0   0   0;    # Na
             0.5 0.5 0.5;  # Cl
             0.5 0   0;    # Cl
             0   0.5 0;    # Cl
             0   0   0.5]  # Cl

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Net ion charges
   Qval = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]
   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha

   @debug "NaCl: Eelec = $E [Ha]"
   # Reference energy from Malik Mamode@DOI 10.1007/s10910-016-0705-9
   @test E/4 ≈ −1.747_564_594_633_182_190_636_212_035
end

@testset "Electrostatic sum: CsCl" begin

   Rd_hat = 2.0 # parameter from Table I

   # CsCl
   a₁ = [1.0, 0.0, 0.0]
   a₂ = [0.0, 1.0, 0.0]
   a₃ = [0.0, 0.0, 1.0]

   r_ions = [0   0.0 0.0;  # Cs; fractional coordinates
             0.5 0.5 0.5]  # Cl

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Net ion charges
   Qval = [1.0, -1.0]
   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha

   @debug "CsCl: Eelec = $E [Ha]"
   @test 0.5E * sqrt(3) ≈ −1.762_674_773_070_99
end


@testset "Electrostatic sum: ZnS" begin

   Rd_hat = 2.0 # parameter from Table I

   # ZnS
   a₁ = [0.0, 0.5, 0.5]
   a₂ = [0.5, 0.0, 0.5]
   a₃ = [0.5, 0.5, 0.0]

   r_ions = [0.00 0.00  0.00;  # Zn; fractional coordinates
             0.25 0.25 0.25]   # S

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Net ion charges
   Qval = [1.0, -1.0]
   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha

   @debug "ZnS: Eelec = $E [Ha]"
   @test 0.25*E * sqrt(3) ≈ −1.638_055_053_388_79
end


@testset "Reproduce NaCl result for IQA solids paper" begin

   Rd_hat = 2.0 # parameter from Table I

   # NaCl
   a₁ = [1.0, 0.0, 0.0] * 10.944
   a₂ = [0.0, 1.0, 0.0] * 10.944
   a₃ = [0.0, 0.0, 1.0] * 10.944

   r_ions = [0   0.5 0.5;  # Na; fractional coordinates
             0.5 0   0.5;  # Na
             0.5 0.5 0;    # Na
             0   0   0;    # Na
             0.5 0.5 0.5;  # Cl
             0.5 0   0;    # Cl
             0   0.5 0;    # Cl
             0   0   0.5]  # Cl

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # QTAIM charges
   Qval = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0] .* 0.88

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @debug "NaCl: Eelec = $E [Ha]"
   @test E ≈ -0.9893 atol=1e-4 # Calculated with ENVIRON
end

@testset "Reproduce CsCl result for IQA solids paper" begin

   Rd_hat = 2.0 # parameter from Table I

   # CsCl
   a₁ = [1.0, 0.0, 0.0] * 7.7199
   a₂ = [0.0, 1.0, 0.0] * 7.7199
   a₃ = [0.0, 0.0, 1.0] * 7.7199

   r_ions = [0   0.0 0.0;  # Cs; fractional coordinates
             0.5 0.5 0.5]  # Cl

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # QTAIM charges
   Qval = [1.0, -1.0] * 0.83
   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha

   @debug "CsCl: Eelec = $E [Ha]"
   @test E ≈ −0.1816 atol=1e-4 # Calculated with ENVIRON
end

@testset "Reproduce LiCl result for IQA solids paper" begin

   Rd_hat = 2.0 # parameter from Table I

   # LiCl
   a₁ = [1.0, 0.0, 0.0] * 9.6933
   a₂ = [0.0, 1.0, 0.0] * 9.6933
   a₃ = [0.0, 0.0, 1.0] * 9.6933

   r_ions = [0   0.5 0.5;  # Li; fractional coordinates
             0.5 0   0.5;  # Li
             0.5 0.5 0;    # Li
             0   0   0;    # Li
             0.5 0.5 0.5;  # Cl
             0.5 0   0;    # Cl
             0   0.5 0;    # Cl
             0   0   0.5]  # Cl

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # QTAIM charges
   Qval = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0] .* 0.90

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
   @debug "LiCl: Eelec = $E [Ha]"
   @test E ≈ -1.1682 atol=1e-4 # Calculated with ENVIRON
end

@testset "Reproduce BN (cubic) result for IQA solids paper" begin

   Rd_hat = 2.0 # parameter from Table I

   # BN cubic
   a₁ = [1.0, 0.0, 0.0] * 6.822
   a₂ = [0.0, 1.0, 0.0] * 6.822
   a₃ = [0.0, 0.0, 1.0] * 6.822

   r_ions = [0   0.5 0.5;  # B; fractional coordinates
             0.5 0   0.5;  # B
             0.5 0.5 0;    # B
             0   0   0;    # B
             0.25 0.25 0.25;  # N
             0.25 0.75 0.75;  # N
             0.75 0.25 0.75;  # N
             0.75 0.75 0.25]  # N

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # QTAIM charges
   Qval = [1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0] .* 2.16

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
	@debug "BN (cubic): Eelec = $E [Ha]"
   @test E ≈ -10.3486 atol=1e-4 # Calculated with ENVIRON
end

@testset "Reproduce BN (hex.) result for IQA solids paper" begin

   Rd_hat = 2.0 # parameter from Table I

   # BN hex
	a₁ = [1.0, -sqrt(3), 0.0] * 4.7325/2
	a₂ = [1.0,  sqrt(3), 0.0] * 4.7325/2
   a₃ = [0.0,  0.0,     1.0] * 12.5834

   r_ions = [0.33333 0.66667 0.25;  # B; fractional coordinates
             0.66667 0.33333 0.75;  # B
             0.33333 0.66667 0.75;  # N
             0.66667 0.33333 0.25]  # N

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # QTAIM charges
   Qval = [1.0, 1.0, -1.0, -1.0] .* 2.214

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat) # Ha
	@debug "BN (hex.): Eelec = $E [Ha]"
   @test E ≈ -5.5359331 atol=1e-6 # Calculated with ENVIRON
end

@testset "Reproduce MgB₂ result for IQA solids paper (3 digits)" begin

   Rd_hat = 2.0 # parameter from Table I

   # MgB₂
	a₁ = [1.0, -sqrt(3), 0.0] * 5.8262/2
	a₂ = [1.0,  sqrt(3), 0.0] * 5.8262/2
   a₃ = [0.0,  0.0,     1.0] * 6.6395

   r_ions = [0.00000 0.00000 0.00;  # Mg; fractional coordinates
             0.33333 0.66667 0.50;  # B
             0.66667 0.33333 0.50]  # B

   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian

   # QTAIM charges
   Qval = [2.0, -1.0, -1.0] .* 0.81

   E = lattice_electrostatic_sum(a₁,a₂,a₃,r_ions,Qval,Rd_hat; shift=1e-3) # Ha
	@debug "MgB₂: Eelec = $E [Ha]"
   @test E ≈ -0.62392 atol=1e-3 # Calculated with ENVIRON
   @test_broken E ≈ -0.62392 atol=1e-4 # Calculated with ENVIRON
   # Warning: the shift is carefully chosen to provide the best estimation,
   # which can only go as close as 1e-3 to the Ewald result
end
