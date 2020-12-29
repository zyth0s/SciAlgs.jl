# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: jl:light,ipynb
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Julia 1.5.3
#     language: julia
#     name: julia-1.5
# ---

# # Hands-on quantum computing for the very curious

# Have you ever heard of quantum teleportation and wondered what is it? Until
# very recently I was like you.  Even if I have received academic formation
# (involving quantum mechanics) I never thought that it was an easy phenomenom
# to grasp - always wondering what the exact protocol is. For instance, once I
# opened a quantum information book at a university public library and I was
# scared enough to not look further. On the other hand, popular explanations of
# the concept always fail to give a satisfactory description because they omit
# the only tool that helps us understand quantum mechanics: math. 
#
# My perspective changed completely when I encountered Nielsen and Matuschak's
# fantastic essays [Quantum Computing for the very
# curious](https://quantum.country/qcvc). Therefore I recommend you to read
# these before continuing. It turns out that
# most of the calculations can be reduced to simple low dimensional matrix
# operations. Furthermore, the whole protocol can be condensed into a single diagram.
#
# In this tutorial we will cover some of the same concepts with the help of a
# computer. All code below is written in Julia.   The linear algebra objects that we need are: the identity matrix
# $I_m$, the adjoint operation $A^\dagger$, the determinant, trace, kroenecker
# (or tensor) product, and a normalization function. The rest is already loaded
# by Julia.
#
#

using LinearAlgebra: I, adjoint, det, tr, kron, normalize!
const ⊗ = kron

# ## Computational quantum basis
#
# Our digital computers are able to reduce all information to zeros and ones,
# that is, a bit may be one of the set $\{0,1\}$. Measurements can be described
# by a probability distribution on the space of possible bits $p: \{0,1\} \to
# [0,1]$. However, quantum bits belong to a two-dimensional complex vector
# space $\mathbb{C}^2$. Measurements give also a probability but draw from a
# larger space $\rho_p: \mathbb{C}^2 \to [0,1]$.
#
# Now you might ask how do we move data stored in bits to qubits. It turns out
# that you can map every classical probability distribution to a quantum
# probability distribution. This distribution is called a density operator,
# $\rho_p$. If we think of each element of the bitset $\{𝟎,𝟏\}$ as a $2\times
# 1$ basis vector, both form an orthonormal basis that spans $\mathbb{C}^2$.
#
# * $|0\rangle$ = [1 0]ᵀ is the classical bit 0
# * $|1\rangle$ = [0 1]ᵀ is the classical bit 1

𝟎 = [1, 0] # ≡ $|0\rangle$; typeset with \bfzero
𝟏 = [0, 1] # ≡ $|1\rangle$; typeset with \bfone

# The density operator is the $2 \times 2$ matrix with the classical
# probabilities arranged along its diagonal, the rest zeros. If we have a
# classical bit 0 (probability $p=1$ to have 𝟎 and $p=0$ to have 𝟏), then the
# density operator is

ρ𝟎 = [1 0;
      0 0]

# Other way to map classical to quantum probabilities is using the outer product

@assert ρ𝟎 ≈ 1𝟎*𝟎' + 0𝟏*𝟏' # 𝟎 has p=1, and 𝟏 has p=0 probability

# The same probability distribution can be arranged as a vector if we had
# chosen the tensor product instead of the outer product. Instead of a $2
# \times 2$ matrix we have a $4 \times 1$ vector. Technically, a vector space
# like $\mathbb{C}^2$ is isomorphic to $\mathbb{C} \otimes \mathbb{C}$.

@assert vec(ρ𝟎) ≈ 1𝟎 ⊗ 𝟎 + 0𝟏 ⊗ 𝟏

# ### A general qubit
#
# A general qubit $|\psi\rangle = \alpha |𝟎\rangle + \beta |𝟏\rangle \in
# \mathbb{C}^2$ is described as a linear combination of the basis vectors and
# it must satisfy the normalization condition $|\alpha|^2 + |\beta|^2 = 1$.

α = 0.6
β = 0.8
@assert α^2 + β^2 ≈ 1 # normalization contraint (here for α, β ∈ ℝ)
ψ = α*𝟎 + β*𝟏
@assert ψ ≈ [α, β]
@assert ψ'*ψ ≈ 1 # ⟨ψ|ψ⟩ = 1 (normalization again)

# ## Quantum gates
#
# Once we have transformed our information to qubits we need to manipulate them
# to achieve our goal. Those operations are realized with gates that are
# analogous to classical circuit gates.
#
# ### Pauli σₓ gate: `X == NOT == ---[ X ]---`

X = [0 1; 1 0]

@assert X*𝟎 ≈ 𝟏
@assert X*𝟏 ≈ 𝟎
# By the linearity of matrix multiplication it follows that the matrix acts the
# same way as the X on all input states, and so they are the same operation.


@assert X*adjoint(X) ≈ I(2) # this is a proof of unitariness

X*ψ

# ### Hadamard gate:  `H == ---[ H ]---`
#
# A Hadamard gate converts a classical bit into a general qubit.

H = [1 1; 1 -1]./√2

@assert H*𝟎 ≈ [1, 1]./√2 # kinda bonding state
@assert H*𝟏 ≈ [1,-1]./√2 # kinda antibonding state
# By the linearity of matrix multiplication it follows that the matrix acts the
# same way as the Hadamard on all input states, and so they are the same
# operation.

@assert H^2 ≈ I(2) # H is idempotent

@assert H*adjoint(H) ≈ I(2) # this is a proof of unitariness
@assert abs(det(H)) ≈ 1 # conservation of particles


# Not every matrix is a gate, for instance

J = [1 1; 1 1]./√2
@assert (J*adjoint(J) ≈ I(2)) == false

# H X ψ == ---[ X ]---[ H ]---
#
# ### Measurement: `---| m )===`
#
# Measurements can be performed by "casting the shadow" of the qubit state at a
# place where we can look, that is, our computational quantum basis. Measuring
# the probability of having the a $|0\rangle$ is the like taking a marginal
# probability from the joint probability.

measure(ψ, verbose=false) = begin
   p𝟎 = tr(ψ*ψ' * 𝟎*𝟎') # ≡ ⟨ψ|0⟩⟨0|ψ⟩ resol. ident
   p𝟏 = tr(ψ*ψ' * 𝟏*𝟏') # ≡ ⟨ψ|1⟩⟨1|ψ⟩ resol. ident
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

#     if m = 0 => input was (|0⟩ + |1⟩)/√2
#     if m = 1 => input was (|0⟩ - |1⟩)/√2
#
# ###  Pauli σy gate:  `Y == ---[ Y ]---`


Y = [0 -im; im 0]
@assert Y*adjoint(Y) ≈ I(2)

# ### Pauli σz gate: `Z == ---[ Z ]---`

Z = [1 0; 0 -1]
@assert Z*adjoint(Z) ≈ I(2)

# ### General rotation gate

θ = π/2

R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

@assert R*adjoint(R) ≈ I(2)

# ## Multi-qubit states
#
# Here again, we can find an analogy with the classical world if we think that
# we are passing from sets to vector spaces. Whereas classical state spaces are
# combined with the cartesian product, such as $\{0,1\} \times \{0,1\}$,
# quantum state spaces are combined with the tensor product, $\mathbb{C}^2 ⊗
# \mathbb{C}^2$. The classical probability $p: \{0,1\} \times \{0,1\} \to
# [0,1]$ on the combined set will be a joint probability (the probability of
# both having 0 in the first and 0 in the second bit, for example). The
# corresponding density operator will be also a joint probability in their
# tensor product space $\rho_p: \mathbb{C}^2 \otimes \mathbb{C}^2 \cong
# \mathbb{C}^4 \to [0,1]$.

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

# ## Multi-qubit gates
# Take the form  G₁ ⊗ G₂ ⊗ ⋯
#
# ### Controlled-NOT gate (CNOT)
#       x ---⋅---
#            |
#       y ---⊕---
#       x is the control qubit
#       y is the target qubit
#       |x, y ⊕ x⟩ for short

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


# Beware, CNOT can change the control qubit!

# |+-⟩
pm = H*𝟎 ⊗ H*𝟏
# |--⟩
mm = H*𝟏 ⊗ H*𝟏

#     |0⟩ ---[ H ]--- |+⟩---⋅--- |-⟩
#                           |
#     |1⟩ ---[ H ]--- |-⟩---⊕--- |-⟩

@assert CNOT*pm ≈ mm

# ## Global phase factor

θ = rand() # any real number

G(θ) = ℯ^(im*θ) * I(2) # global phase factor ℯ^(iθ)

if typeof(θ) <: Real
   @assert G(θ)*adjoint(G(θ)) ≈ I(2)
   @info "A matrix changing the global phase factor is unitary."
end

@info "Changing the phase does not affect the measurement."
@assert G(rand())*𝟎 |> measure ≈ measure(𝟎)

# ## Other gates

S = [ 1 0; 0 im]
T = [ 1 0; 0 ℯ^(im*π/4)]
#Y = [ 0 -im; im 0]
#Z = [ 1 0; 0 -1]

# ## Quantum teleportation

# Alice and Bob can share beforehand a special two-qubit
#
#     |0⟩ ---[ H ]---⋅--- 
#                    |   
#                    |    (|00⟩ + |11⟩)/√2
#                    |
#     |0⟩ -----------⊕---
#     
#


ebit = CNOT*H₁*𝟎𝟎 # entangled bit -> shared

# Then, they can separate.
#
# The full protocol can be condensed in the following diagram
#
#              teleported state:  |ψ⟩  ------⋅---[ H ]---[ z )===
#                                            |
#                                            |
#                                            |
#            |0⟩ ---[ H ]---⋅----------------⊕-----------[ x )===
#                           |  |00⟩ + |11⟩
#                           | ------------
#                           |      √2
#            |0⟩ -----------⊕------------------------------------[ Xˣ ]---[ Zᶻ ]--- |ψ⟩
#            
# Given the moment, Alice wants to send Bob some important piece of information
# encoded as a qubit.

# Any state ψ we want to teleport
α = rand(Complex{Float64})
β = sqrt(1 - α*conj(α)) # |α|² + |β|² = 1
@assert α*conj(α) + β*conj(β) ≈ 1 "State not properly normalized. Try with other (α,β)"
ψ = α*𝟎 + β*𝟏 # ∈ ℂ² ≝ ℂ ⊗ ℂ ; ρψ = ψ*ψ' ∈ ℂ ⊗ ℂ
_ψ = ψ # we can do this only in a classic circuit (debugging purposes)

# She applies a conditional-NOT to the entangled qubit based on the state to be
# teleported, and applies a Hadamard matrix to the state to be teleported.

s = ψ ⊗ ebit;

gate1 = CNOT ⊗ I(2)
gate2 = H ⊗ I(4);

ψ = gate2*gate1*s;

# Alice measures first two bits, posibilities are 00, 01, 10, and 11

P𝟎𝟎 = 𝟎𝟎*𝟎𝟎' ⊗ I(2) # projections
P𝟎𝟏 = 𝟎𝟏*𝟎𝟏' ⊗ I(2)
P𝟏𝟎 = 𝟏𝟎*𝟏𝟎' ⊗ I(2)
P𝟏𝟏 = 𝟏𝟏*𝟏𝟏' ⊗ I(2);

ρψ = ψ*ψ' # density operator
p𝟎𝟎 = tr(ρψ * P𝟎𝟎) |> real # probabilities
p𝟎𝟏 = tr(ρψ * P𝟎𝟏) |> real
p𝟏𝟎 = tr(ρψ * P𝟏𝟎) |> real
p𝟏𝟏 = tr(ρψ * P𝟏𝟏) |> real;

@info "  The probability of |𝟎𝟎⟩ is $p𝟎𝟎"
@info "  The probability of |𝟎𝟏⟩ is $p𝟎𝟏"
@info "  The probability of |𝟏𝟎⟩ is $p𝟏𝟎"
@info "  The probability of |𝟏𝟏⟩ is $p𝟏𝟏"

# All outcomes have the same chance.

icollapsed = argmax([p𝟎𝟎, p𝟎𝟏, p𝟏𝟎, p𝟏𝟏])
icollapsed = rand(1:4) # to avoid taking always the first
Pcollapsed = [P𝟎𝟎, P𝟎𝟏, P𝟏𝟎, P𝟏𝟏][icollapsed];

x = (icollapsed == 2 || icollapsed == 4) |> Int
z = (icollapsed == 3 || icollapsed == 4) |> Int
@info "  Alice measured x = $x and z = $z"
@info "  Alice qubits collapsed to $(["|𝟎𝟎⟩", "|𝟎𝟏⟩", "|𝟏𝟎⟩", "|𝟏𝟏⟩"][icollapsed])"

# The resulting two classical bits are shared to Bob, who uses them to process
# his entangled bit.

range = 2icollapsed-1:2icollapsed
ψ = normalize!(Pcollapsed*ψ)[range] # state after measurement of 𝟎𝟎
ψ = Z^z * X^x * ψ # Bob uses Alice classical bits x and z
@assert _ψ ≈ ψ "Teleported state has been corrupted"
@info "  Teleported |ψ⟩ = ($(ψ[1])) |𝟎⟩ + ($(ψ[2])) |𝟏⟩ !!"


# This is awesome! We have teleported a vast amount information just moving two
# bits. It is like seeing a sportman moving a truck with his beard! 

# ## Toffoli gate CCNOT

#CCNOT = cat(I(6), [0 1; 1 0], dims=(1,2)) |> Matrix
CCNOT = cat(I(4), CNOT, dims=(1,2)) |> Matrix

Rx(θ) = cos(θ/2)*I(2) - im * sin(θ/2)*X
Ry(θ) = cos(θ/2)*I(2) - im * sin(θ/2)*Y
Rz(θ) = cos(θ/2)*I(2) - im * sin(θ/2)*Z

