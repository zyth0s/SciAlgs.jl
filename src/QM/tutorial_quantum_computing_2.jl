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

# # Encoding Ising-like Hamiltonians as quantum circuits
#
# A more lengthy explanation is given in
# https://www.mustythoughts.com/variational-quantum-eigensolver-explained
using LinearAlgebra: I, tr

# Computational basis
𝟎 = [1,0]
𝟏 = [0,1]

X = [0 1; 1 0]
Z = [1 0; 0 -1] # we measure in Z basis

H = [1 1; 1 -1]/√2
p = H * 𝟎 # = |+⟩
m = H * 𝟏 # = |-⟩

# ## Measurements in fixed Z-basis: rotate to measure other Pauli matrices
#
# If we measure in the Z basis its not possible to distinguish between |+⟩ and
# |-⟩. 

@info "⟨+|Z|+⟩ ≈ ⟨-|Z|-⟩", p'*Z*p ≈ m'* Z*m

# However, if we rotate the basis, then it is possible.

Ry(θ) = [cos(θ/2) -sin(θ/2); sin(θ/2) cos(θ/2)]

@assert Ry(π/2)*p ≈ 𝟏
@assert Ry(π/2)*m ≈ 𝟎

p_rot = Ry(-π/2)*p
m_rot = Ry(-π/2)*m
@info "⟨+|Z|Ry(-π/2)|+⟩ ≈ ⟨-|Z|Ry(-π/2)|-⟩", p_rot'*Z*p_rot ≈ m_rot'* Z*m_rot
@info "Outcome of ⟨+|Z|Ry(-π/2)|+⟩ = ", p_rot'*Z*p_rot
@info "Outcome of ⟨-|Z|Ry(-π/2)|-⟩ = ", m_rot'* Z*m_rot

# The solution is to apply
# * $R_y(−\pi/2)$      if Hamiltonian has $\hat{X}$
# * $R_x(\pi/2)$       if Hamiltonian has $\hat{Y}$
# * $I$                if Hamiltonian has $\hat{Z}$
#
# ## Divide the Hamiltonian in single Pauli terms

𝓗₁ = 2Z
𝓗₂ = X
𝓗₃ = I(2)
𝓗 = 𝓗₁ + 𝓗₂ + 𝓗₃
@assert 𝓗 ≈ [3 1; 1 -1]

# ## Initial state
#
# We use a parametrized initial wavefunction for reasons that will be explained
# later but you can already guess that we will do a variational search.

# Ansatz will be Ry(θ)*𝟎
ansatz(θ) = Ry(θ)*𝟎

# ## A circuit for each term

for θ in [0,π]

   @info "θ = $θ"
   # Circuit for 𝓗₁: ψ ---[ h₁ )===

   ψ = ansatz(θ)

   h₁ = 2ψ' * Z * ψ
   @info "  E₁ = $h₁"

   # Circuit for 𝓗₂: ψ ---[ Ry ]---[ h₂ )===

   ψ = Ry(-π/2)*ansatz(θ)

   h₂ = ψ' * Z * ψ
   @info "  E₂ = $h₂"

   # Circuit for 𝓗₃:

   h₃ = 1
   @info "  E₃ = $h₃"

   @info "  ⟨ψ($θ)|𝓗|ψ($θ)⟩ = $(h₁ + h₂ + h₃)"
end
