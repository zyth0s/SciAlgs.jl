# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.1
#   kernelspec:
#     display_name: Julia 1.0.2
#     language: julia
#     name: julia-1.0
# ---

# ## Spectral decomposition with Julia

# We discuss the diagonalization of a square matrix to obtain
# its eigenvalues and eigenvectors.

# Note about notation: 
# * vectors have a bar above, like $v̄$.
# * matrices have a bar below, like $A̲$.
# * scalars have no bar.
# * ' means (i) the complex cojugate for scalars, or (ii) the (conjugate) transpose for vectors/matrices.

using LinearAlgebra: eigen, Diagonal, I, det, rank

# Lets take a $4 \times 4$ matrix

A = rand(4,4)


# that is symmetric (or hermitic)

A = (A + A')/2.

# A $4 \times 4$ matrix has 4 eigenvalues $(e_1,e_2,e_3,e_4)$
# and 4 respective eigenvectors $(v̄_1,v̄_2,v̄_3,v̄_4)$

e, U = eigen(A) # solves $A̲v̄ = ev̄$

# `eigen` returns the eigenvalues in a vector `e` and
# the eigenvectors as columns of a matrix `U`.
# E.g. the first eigenvalue is `e1 = e[1]` and the first 
# eigenvector is `v̄1 = U̲[:,1]`

# We can check that the secular equation $|A̲ - eI̲| = 0$
# is satisfied (for each eigenvalue)

for i in 1:4
   isapprox( det(A - e[i]*I), 0.00, atol=1e-4) || error("Secular equation $i not satisfied")
end

# The corresponding eigenvectors satisfy $(A̲ - eI̲)v̄ = 0̄$

for i in 1:4
   isapprox( (A - e[i]*I)*U[:,i], fill(0.00,4), atol=1e-4) || error("a) Eigenvalue equation $i not satisfied")
end

# or equivalently $A̲v̄ = ev̄$

for i in 1:4
   isapprox( A*U[:,i], e[i]*U[:,i], atol=1e-4) || error("b) Eigenvalue equation $i not satisfied")
end

# - If $A̲$ is symmetric, thus $U̲$ is orthonormal and $U̲' = U̲^{-1}$
# - If $A̲$ is hermitic, thus $U̲$ is unitary and $U̲' = U̲^{-1}$
#
# In both cases it can be described as a similarity transformation
#  $$ A̲ = U̲ D̲ U̲^{-1} = U̲ D̲ U̲' $$
# where $D̲$ is a diagonal matrix with the eigenvalues as diagonal elements.
# Therefore, $U̲$ is the matrix that diagonalizes $A̲$.
# $A̲$ and $D̲$ are said to be similar matrices.

A ≈ U * Diagonal(e) * U' || error("Not a similarity transformation")

# A similarity transformation is also a conformal mapping. A matrix $W̲$
# encoding a crystallographic symmetry operation is an isometry if
# $$ G̲ ≈ W̲ G̲ W̲' $$
# where $G̲$ is the metric tensor.

# ---

# # Eigenvalues and eigenvectors: Generalized and determinant free

# Jeff Suzuki (2020)
# Eigenvalues and Eigenvectors: Generalized and Determinant Free
# Mathematics Magazine, 93:3, 200-212, DOI: 10.1080/0025570X.2020.1736874

using RowEchelon
#using Polynomials (unable to precompile error)

# <div class="alert alert-block alert-danger">
# <b>Error:</b> Precompilation of Polynomials failed.
# </div>

# The matrix to be diagonalized is:

N = 3
#A = rand(N,N)
A = [3 4 8; 8 7 16; -4 -4 -9]

# There are exactly $N$ eigenvalues and independent eigenvectors

# ## Eigenvalues

# 1. Choose a random non zero vector $u$

#u = rand(N)
u = [1, 0, 0]
@assert u != zeros(N)
u


# 2. Find $k$ such that all elements of
#    $V_k = (u, A u, A^2 u, A^3 u, ..., A^k u)$ are linearly independent.

Vk = zeros(N,N+1)
Vk[:,1] = u
for k in 1:N
  Vk[:,k+1] = A^k*u
end
Vk

# Gaussian elimination of V_k to row echelon format:

Vk_ref = rref(Vk)

k = rank(Vk)
println("The minimum polynomial is of degree $k.")
println("        A^$k u + ",join(["x_$i A^$i u + " for i in k:-1:2]),"x_1 u = 0")
println("f(A) =  A^$k   + ",join(["x_$i A^$i u + " for i in k:-1:2]),"x_1 I = 0")

# <div class="alert alert-block alert-warning">
# <b>Warning:</b> The minimal potential is not found automatically. Needs human intervention.
# </div>

# By visual inspection we see that 3u + 2(A u) = A^2 u.
# x_1 = -3, x_2 = -2.
# The coefficients of the minimal polynomial are

pmin = [-3, -2, 1]


# Knowing that the minimal polynomial accepts this form

println("f(A) = ",join(["(A-λ_$i I)" for i in 1:k]))

# eigvalues = roots(Polynomial(pmin)) # reads e.g. P + P^1 + P^2 + ... = 0

# ## Eigenvectors

# rref(A - Diagonal(fill(eigvalues[1],N)))
# rref(A - Diagonal(fill(eigvalues[2],N)))

