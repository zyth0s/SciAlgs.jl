
# 10.1119/10.0002644

using LinearAlgebra

function build_c_F(ϑ; Z=1)
   # Orbital energies [Ha]. Table I.
   ϵ = [-Z^2/2, -Z^2/8]
   # Coulomb integrals [Ha]. Table I.
   I₁₁₁₁ = 5/8 * Z
   I₁₁₁₂ = 2^12 * √2/27/7^4 * Z
   I₁₁₂₂ = 16/9^3 * Z
   I₁₂₁₂ = 17/3^4 * Z
   I₁₂₂₂ = 2^9 * √2/27/5^5 * Z
   I₂₂₂₂ = 77/2^9 * Z

   c = [cos(ϑ), sin(ϑ)] # eq 16 & 17
   # Fock matrix elements
   F = zeros(2,2)
   F[1,1] = ϵ[1] +   I₁₁₁₁ * c[1]^2 + 2I₁₁₁₂ * c[1] * c[2] + I₁₂₁₂ * c[2]^2 # eq 23
   F[1,2] = F[2,1] = I₁₁₁₂ * c[1]^2 + 2I₁₁₂₂ * c[1] * c[2] + I₁₂₂₂ * c[2]^2 # eq 24 & 25
   F[2,2] = ϵ[2] +   I₂₂₂₂ * c[2]^2 + 2I₁₂₂₂ * c[1] * c[2] + I₁₂₁₂ * c[1]^2 # eq 26

   # Calculate two-electron energy eq 13 & 29
   Egs = tr(Diagonal(ϵ)*c*c') + tr(F*c*c')
   @assert Egs ≈ 2tr(Diagonal(ϵ)*c*c') +
                  I₁₁₁₁ * c[1]^4 +
                  I₂₂₂₂ * c[2]^4 +
                 4I₁₁₁₂ * c[1]^3 * c[2] +
                 4I₁₂₂₂ * c[2]^3 * c[1] +
                 2(2I₁₁₂₂ + I₁₂₁₂) * c[1]^2 * c[2]^2

   c, F, ϵ, Egs
end

# Initialization of the HF SCF procedure
let Z = 2,         # Atomic number
    max_iter = 20, # Maximum number of SCF iteractions
    ϑ = 0,         # Initial polar angle
    Egs = 0.0      # Initial total [ground state] energy

   # Iterative SCF loop
   for i in 1:max_iter
      c, F, ϵ, Egs = build_c_F(ϑ; Z)

      # Calculate lower eigenvalue of the Fock matrix ([canonical] orbital energy)
      ε = eigvals(F) |> first
      # check against standard quadratic formula
      @assert ε ≈ 0.5(F[1,1] + F[2,2]) - sqrt(0.25(F[1,1] - F[2,2])^2 + F[1,2]*F[2,1])

      @info "$i ϑ = $ϑ"
      @info "   Total   energy = $Egs [Ha]"
      @info "   Orbital energy = $ε [Ha]"

      # Polar angle
      ϑ = atan((ε-F[1,1])/F[1,2])

   end
   Z == 2 && @assert isapprox(Egs, -2.82, atol=1e-2) "$Egs"
   Z == 2 && @assert Egs ≈ -2.82363522301439
end

# TODO: Egs vs. ϑ
# INIT: 4πr²χ₁ vs. r  & 4πr²χ₂ vs. r  &  4πr²φHF vs. r
# TODO: ϵ(HF) vs. Z  &  ϑ(HF) vs. Z
# TODO: Egs(HF) vs. i  &  ϑ(HF) vs. i

let a₀ = 1 # 5.29177210903e-11 #[m]
   χ₁(r,Z) = inv(√π) * (Z/a₀)^(3/2)                   * exp(-Z*r/a₀)
   χ₂(r,Z) = inv(√π) * (Z/a₀)^(3/2) * (1 - Z*r/(2a₀)) * exp(-Z*r/(2a₀))
   using PyPlot
   r = 0:0.01:5
   Z = 2
   plot(r,4π*r.^2 .* χ₁.(r,Z).^2, "--", label="4πr² χ₁(r)²")
   plot(r,4π*r.^2 .* χ₂.(r,Z).^2, "--", label="4πr² χ₂(r)²")
   φHF = (cos(-0.2)* χ₁.(r,Z) + sin(-0.2)*χ₂.(r,Z))
   plot(r,4π*r.^2 .* φHF.^2, "-", label="4πr² φHF(r)²")
   ylim(0,4)
   xlim(0,5)
   legend()
end
