# Quantum computing for the very curious

using LinearAlgebra: I, adjoint, det, tr, kron, normalize!
const ⊗ = kron
@info "* Computational quantum basis"
@info "  |0⟩ = [1 0]ᵀ is the classical bit 0"
@info "  |1⟩ = [0 1]ᵀ is the classical bit 1"
𝟎 = [1, 0] # \bfzero
𝟏 = [0, 1] # \bfone

@info "* A qubit"
@info "  |ψ⟩ = α |𝟎⟩ + β |𝟏⟩"
α = 0.6
β = 0.8
@info "  Normalization constraint: |α|² + |β|² = 1"
@assert α^2 + β^2 ≈ 1 # normalization contraint (for α, β ∈ ℝ)
ψ = α*𝟎 + β*𝟏
@assert ψ ≈ [α, β]

@info "* Quantum gates"
@info "  - NOT == X == σₓ == ---[ X ]---"

X = [0 1; 1 0]

@assert X*𝟎 ≈ 𝟏
@assert X*𝟏 ≈ 𝟎
# By the linearity of matrix multiplication it follows that the matrix acts the
# same way as the X on all input states, and so they are the same
# operation.


@assert X*adjoint(X) ≈ I(2) # this is a proof of unitariness

X*ψ

@info "  - Hadamard gate == H == ---[ H ]---"

H = [1 1; 1 -1]./√2

@assert H*𝟎 ≈ [1, 1]./√2 # kinda bonding state
@assert H*𝟏 ≈ [1,-1]./√2 # kinda antibonding state
# By the linearity of matrix multiplication it follows that the matrix acts the
# same way as the Hadamard on all input states, and so they are the same
# operation.

@assert H^2 ≈ I(2)

@assert H*adjoint(H) ≈ I(2) # this is a proof of unitariness
@assert abs(det(H)) ≈ 1 # conservation of particles


#J = [1 1; 1 1]./√2
#@assert abs(det(J)) ≈ 1 # check this is a proof of non-unitariness

# H X ψ == ---[ X ]---[ H ]---

@info "  Measurement == ---| m )==="
measure(ψ, verbose=false) = begin
   #println("m = 0 with probability $(ψ[1]^2)")
   #println("m = 1 with probability $(ψ[2]^2)")
   #p𝟎 = ψ[1]*conj(ψ[1])
   #p𝟏 = ψ[2]*conj(ψ[2])
   p𝟎 = tr(ψ*ψ' * 𝟎*𝟎')
   p𝟏 = tr(ψ*ψ' * 𝟏*𝟏')
   if verbose
      println("m = 0 with probability $p𝟎")
      println("m = 1 with probability $p𝟏")
   end
   [p𝟎, p𝟏]
end

# _ψ2 ---[ H ]---[ m )===
_ψ = rand([𝟎,𝟏])
_ψ2 = H*_ψ # Input state is either (|0⟩ + |1⟩)/√2 or (|0⟩ - |1⟩)/√2
H*_ψ2 |> measure

# if m = 0 => input was (|0⟩ + |1⟩)/√2
# if m = 1 => input was (|0⟩ - |1⟩)/√2


@info "  - Y == σy == ---[ Y ]---"
Y = [0 -im; im 0]
@assert Y*adjoint(Y) ≈ I(2)

@info "  - Z == σz == ---[ Z ]---"
Z = [1 0; 0 -1]
@assert Z*adjoint(Z) ≈ I(2)

@info "  - General rotation"

θ = π/2

R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

@assert R*adjoint(R) ≈ I(2)

@info "* Multi-qubit states"
@info "  ψ₁ ⊗ ψ₂ ⊗ ⋯"

𝟎𝟎 = 𝟎 ⊗ 𝟎
𝟎𝟏 = 𝟎 ⊗ 𝟏
𝟏𝟎 = 𝟏 ⊗ 𝟎
𝟏𝟏 = 𝟏 ⊗ 𝟏

γ = 0.8
δ = 0.6
ϕ = γ*𝟎 + δ*𝟏

# More generally, if we have single-qubit states ψ and ϕ, then the combined
# state when the two qubits are put together is just:

ξ = ψ ⊗ ϕ
@assert ξ ≈ [ψ[1]*ϕ[1], ψ[1]*ϕ[2], ψ[2]*ϕ[1], ψ[2]*ϕ[2]]

@info "* Multi-qubit gates"
@info "  G₁ ⊗ G₂ ⊗ ⋯"
@info "  - Controlled-NOT == CNOT"
@info "      x ---⋅---"
@info "           |"
@info "      y ---⊕---"
@info "      x is the control qubit"
@info "      y is the target qubit"
@info "      |x, y ⊕ x⟩ for short"

CNOT = [1 0 0 0;
        0 1 0 0;
        0 0 0 1;
        0 0 1 0]
# Also
@assert CNOT ≈ cat(I(2), X, dims=(1,2)) #|> Matrix
@assert CNOT ≈ 𝟎*𝟎' ⊗ I(2) + 𝟏*𝟏' ⊗ X


@assert CNOT*𝟎𝟎 ≈ 𝟎𝟎
@assert CNOT*𝟎𝟏 ≈ 𝟎𝟏
@assert CNOT*𝟏𝟎 ≈ 𝟏𝟏
@assert CNOT*𝟏𝟏 ≈ 𝟏𝟎


# Apply H to first qubit in a 2d space
H₁ = H ⊗ I(2)
# Apply H to second qubit in a 2d space
H₂ = I(2) ⊗ H

CNOT*H₁*𝟎𝟎


@info "    CNOT can change the control qubit!"
# |+-⟩
pm = H*𝟎 ⊗ H*𝟏
# |--⟩
mm = H*𝟏 ⊗ H*𝟏

# |0⟩ ---[ H ]--- |+⟩---⋅--- |-⟩
#                       |
# |1⟩ ---[ H ]--- |-⟩---⊕--- |-⟩

@assert CNOT*pm ≈ mm



# Global phase factor

θ = rand() # any real number

G(θ) = ℯ^(im*θ) * I(2) # global phase factor ℯ^(iθ)

if typeof(θ) <: Real
   @assert G(θ)*adjoint(G(θ)) ≈ I(2)
   @info "A matrix changing the global phase factor is unitary."
end

@info "Changing the phase does not affect the measurement."
@assert G(rand())*𝟎 |> measure ≈ measure(𝟎)

# Other gates

S = [ 1 0; 0 im]
T = [ 1 0; 0 ℯ^(im*π/4)]
#Y = [ 0 -im; im 0]
#Z = [ 1 0; 0 -1]

@info "* Quantum teleportation"

# Special two-qubit shared between Alice and Bob
# |0⟩ ---[ H ]---⋅--- 
#                |   
#                |    (|00⟩ + |11⟩)/√2
#                |
# |0⟩ -----------⊕---
ebit = CNOT*H₁*𝟎𝟎 # entangled bit -> shared


@info """

             teleported state:  |ψ⟩  ------⋅---[ H ]---[ z )===
                                           |
                                           |
                                           |
           |0⟩ ---[ H ]---⋅----------------⊕-----------[ x )===
                          |  |00⟩ + |11⟩
                          | ------------
                          |      √2
           |0⟩ -----------⊕------------------------------------[ Xˣ ]---[ Zᶻ ]--- |ψ⟩
"""

# Any state ψ we want to teleport
α = rand(Complex{Float64})
β = sqrt(1 - α*conj(α)) # |α|² + |β|² = 1
@assert α*conj(α) + β*conj(β) ≈ 1 "State not properly normalized. Try with other (α,β)"
ψ = α*𝟎 + β*𝟏 # ∈ ℂ² ≝ ℂ ⊗ ℂ ; ρψ = ψ*ψ' ∈ ℂ ⊗ ℂ
_ψ = ψ # we can do this only in a classic circuit (debugging purposes)

s = ψ ⊗ ebit

gate1 = CNOT ⊗ I(2)
gate2 = H ⊗ I(4)

ψ = gate2*gate1*s

# Alice measures first two bits, posibilities: 00, 01, 10, and 11

P𝟎𝟎 = 𝟎𝟎*𝟎𝟎' ⊗ I(2) # projections
P𝟎𝟏 = 𝟎𝟏*𝟎𝟏' ⊗ I(2)
P𝟏𝟎 = 𝟏𝟎*𝟏𝟎' ⊗ I(2)
P𝟏𝟏 = 𝟏𝟏*𝟏𝟏' ⊗ I(2)

ρψ = ψ*ψ' # density operator
p𝟎𝟎 = tr(ρψ * P𝟎𝟎) |> real # probabilities
p𝟎𝟏 = tr(ρψ * P𝟎𝟏) |> real
p𝟏𝟎 = tr(ρψ * P𝟏𝟎) |> real
p𝟏𝟏 = tr(ρψ * P𝟏𝟏) |> real

@info "  The probability of |𝟎𝟎⟩ is $p𝟎𝟎"
@info "  The probability of |𝟎𝟏⟩ is $p𝟎𝟏"
@info "  The probability of |𝟏𝟎⟩ is $p𝟏𝟎"
@info "  The probability of |𝟏𝟏⟩ is $p𝟏𝟏"

icollapsed = argmax([p𝟎𝟎, p𝟎𝟏, p𝟏𝟎, p𝟏𝟏])
icollapsed = rand(1:4) # to avoid taking always the first
Pcollapsed = [P𝟎𝟎, P𝟎𝟏, P𝟏𝟎, P𝟏𝟏][icollapsed]

x = (icollapsed == 2 || icollapsed == 4) |> Int
z = (icollapsed == 3 || icollapsed == 4) |> Int
@info "  Alice measured x = $x and z = $z"
@info "  Alice qubits collapsed to $(["|𝟎𝟎⟩", "|𝟎𝟏⟩", "|𝟏𝟎⟩", "|𝟏𝟏⟩"][icollapsed])"

range = 2icollapsed-1:2icollapsed
ψ = normalize!(Pcollapsed*ψ)[range] # state after measurement of 𝟎𝟎
ψ = Z^z * X^x * ψ # Bob uses Alice classical bits x and z
@assert _ψ ≈ ψ "Teleported state has been corrupted"
@info "  Teleported |ψ⟩ = ($(ψ[1])) |𝟎⟩ + ($(ψ[2])) |𝟏⟩ !!"


@info "Toffoli gate CCNOT"

#CCNOT = cat(I(6), [0 1; 1 0], dims=(1,2)) |> Matrix
CCNOT = cat(I(2), CNOT, dims=(1,2)) |> Matrix
