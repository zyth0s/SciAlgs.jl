
using LinearAlgebra: det, norm
#var"'ᵀ"(a) = transpose(a) # ≥ 1.6 (Julia)

erf( x) = @ccall erf( x::Float64)::Float64
erfc(x) = @ccall erfc(x::Float64)::Float64

@doc """
    lattice_electrostatic_sum(a₁, a₂, a₃, r_ions, ion_charges, Rd_hat=2.0)

Compute the electrostatic sum of an extended distribution of point charges with
Pickard's algorithm (doi:10.1103/PhysRevMaterials.2.013806). The unit cell is
defined by the parallelepiped with vectors `a₁`, `a₂`, `a₃`. Ion coordinates
are given as rows of `r_ions`, and their charges are listed in `ion_charges`.
Parameter `Rd_hat` is a factor for the damping (and real-space) cutoff radius.
It is 2.0 by default.
"""
function lattice_electrostatic_sum(a₁, a₂, a₃, r_ions, ion_charges, Rd_hat=2.0; shift=abs(minimum(ion_charges)) + 1)

   @assert size(a₁) == size(a₂) == size(a₃) == (3,)
   nions = size(r_ions, 1)
   @assert size(r_ions) == (nions,3)
   @assert size(ion_charges) == (nions,)

   # No charge => no energy
   sum(abs.(ion_charges)) ≈ 0 && return 0.0

   # Note that Q must be ≠ 0 for this method to work. To calculate
   # total elec. energy we shift the charge and do the calculation
   # separately for the background
   # initial energy
   E = if sum(ion_charges) ≈ 0
          # Shift charges with a background to have net charge ≠ 0
          #shift = abs(minimum(ion_charges)) + 1
          #shift = 0.001
          #any(q -> q ≈ shift, ion_charges) && error("Shift charge neutralizes the charge of an ion.")
          ion_charges .+= shift
          back_charges = shift*ones(size(ion_charges))
          -lattice_electrostatic_sum(a₁, a₂, a₃, r_ions, back_charges, Rd_hat)
       else
          0
   end

   # Lattice vectors; by the way, the metric tensor is Lᵀ⋅L
   L = hcat(a₁,a₂,a₃)
   # cell volume
   vol = abs(det(L))
   # average density
   ρ = sum(ion_charges) / vol

   # reciprocal lattice vectors
   # [𝐛₁𝐛₂𝐛₃]ᵀ = 2π [𝐚₁𝐚₂𝐚₃]⁻¹
   B = inv(L) # * 2π

   # interplanar distances
   # dₕₖₗ = 2π / |𝐠ₕₖₗ| where 𝐠ₕₖₗ = h 𝐛₁ + k 𝐛₂ + l 𝐛₃
   d₁₀₀ = inv(norm(B[1,:])) # * 2π (cancels)
   d₀₁₀ = inv(norm(B[2,:])) # * 2π (cancels)
   d₀₀₁ = inv(norm(B[3,:])) # * 2π (cancels)

   # Length scale (largest perpendicular distance between the faces of the primitive cell)
   hₘₐₓ = max(d₁₀₀,d₀₁₀,d₀₀₁) + eps(d₀₀₁)
   @debug "" hₘₐₓ
   Rd = Rd_hat * hₘₐₓ
   Rc = 3Rd_hat^2 * hₘₐₓ # optimal real-space cutoff for convergence (caption Fig 2)

   # number of cells included in each direction [hkl]
   ncells₁₀₀ = ceil(Int, Rc / d₁₀₀)
   ncells₀₁₀ = ceil(Int, Rc / d₀₁₀)
   ncells₀₀₁ = ceil(Int, Rc / d₀₀₁)

   # Run over ions in (central) unit cell
   for i in 1:nions

      # initial energy of ion `i`
      Ei = 0.0
      Qi = ion_charges[i] # we will skip the i==j terms

      # run over cells in supercell
      for n₁ in -ncells₁₀₀:ncells₁₀₀,
          n₂ in -ncells₀₁₀:ncells₀₁₀,
          n₃ in -ncells₀₀₁:ncells₀₀₁

         # Current (neighbor) cell coordinates
         Rcell = n₁*a₁ + n₂*a₂ + n₃*a₃

         # Loop over ions in current (neighbor) cell
         for j in 1:nions
            i == j && n₁ == n₂ == n₃ == 0 && continue

            # distance between ions
            rᵢⱼ = norm(r_ions[i,:] - (Rcell + r_ions[j,:]))

            rᵢⱼ > Rc && continue

            Ei += ion_charges[j] * erfc(rᵢⱼ/Rd) / rᵢⱼ # [1st term in eq 13] * 2 / ion_chargesᵢ
            Qi += ion_charges[j]
         end
      end

      Ei *= 0.5ion_charges[i] # finish [eq 13]

      # ΔEᵢ = ΔE(sphere) + ΔE(damp) + ΔE(self) [eq 14]
      Ra = cbrt(3Qi/(4π*ρ)) # [eq 10]
      Ei -= π*ion_charges[i] * ρ * Ra^2                         # ΔE(sphere) [eq 15]
      Ei += π*ion_charges[i] * ρ * (Ra^2 - Rd^2/2) * erf(Ra/Rd) +
           √π*ion_charges[i] * ρ * Ra * Rd * exp(-Ra^2/(Rd^2))  # ΔE(damp)   [eq 16]
      Ei -= inv(√π*Rd) * ion_charges[i]^2                       # ΔE(self)   [eq 18]
      E += Ei
   end
   E
end


# Test: Al; ICSD:43423

function test1()
   Rd_hat = 2.0 # value from Table I

   a₁ = [5.41141973394663 , 0.0              , 0.0]
   a₂ = [2.70570986697332 , 4.68642696013821 , 0.0]
   a₃ = [2.70570986697332 , 1.56214232004608 , 4.41840571073226]

   r_ions = zeros(1,3)
   # Valence charge
   Z = [3.0]

   E = lattice_electrostatic_sum(a₁,a₂,a₃, r_ions, Z, Rd_hat)
   @debug "Al: Eᴺᴺ = $E [Ha]"
   @assert E ≈ -2.695954572 # from Table I
end

# Test: Si; ICSD:51688

function test2()
   Rd_hat = 2.0 # value from Table I

   a₁ = [7.25654832321381, 0.00000000000000, 0.00000000000000]
   a₂ = [3.62827416160690, 6.28435519169252, 0.00000000000000]
   a₃ = [3.62827416160690, 2.09478506389751, 5.92494689524090]

   r_ions = [0.00  0.00  0.00; # Si ; fractional coords
             0.25  0.25  0.25] # Si
   #r_ions = (hcat(a₁,a₂,a₃)*r_ions'ᵀ)'ᵀ # cartesian
   r_ions = (hcat(a₁,a₂,a₃)*r_ions')' # cartesian
   # Valence charges
   Z = [4.0, 4.0]

   E = lattice_electrostatic_sum(a₁,a₂,a₃, r_ions, Z, Rd_hat)
   @debug "Si: Eᴺᴺ = $E [Ha]"
   @assert E ≈ -8.398574646 # from Table I
end

# Test: SiO₂; ICSD:29122

function test3()
   Rd_hat = 2.0 # value from Table I

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
   Z = vcat(repeat([6.0], 6), repeat([4.0], 3))

   E = lattice_electrostatic_sum(a₁,a₂,a₃, r_ions, Z, Rd_hat)
   @debug "SiO₂: Eᴺᴺ = $E [Ha]"
   @assert E ≈ -69.488098659 # from Table I
end

# Test: Al₂SiO₅; ICSD:24275

function test4()
   Rd_hat = 2.0 # value from Table I

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
   Z = 6.0ones(size(r_ions)[1]) # O atoms (initial)
   Z[[9,10,11,12,13,15,17,19]] .= 3.0 # Al atoms
   Z[[21,24,27,30]]            .= 4.0 # Si atoms

   E = lattice_electrostatic_sum(a₁,a₂,a₃, r_ions, Z, Rd_hat)
   @debug "Al₂SiO₅: Eᴺᴺ = $E [Ha]"
   @assert E ≈ -244.055008450 # from Table I
end
test1()
test2()
test3()
test4()
