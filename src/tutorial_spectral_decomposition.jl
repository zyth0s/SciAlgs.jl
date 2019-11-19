
# Script to explain the spectral decomposition,
# also called eigenvalue equation
#
# Note about notation: - vector have a bar above, like v̂
#                      - matrices have a bar below, like A̲
#                      - scalars have no bar
#                      - ' means the (conjugate) transpose


using LinearAlgebra


# Lets take a 4x4 matrix
A = rand(4,4)
# that is symmetric/hermitic
A = (A + A')/2.

# A 4x4 matrix has 4 eigenvalues (e1,e2,e3,e4) 
# and 4 respective eigenvectors (v̄1,v̄2,v̄3,v̄4)
e, U = eigen(A) # solves A̲v̄ = ev̄

# "eigen" returns the eigenvalues in a vector e and
# the eigenvectors as columns of a matrix U̲.
# E.g. the first eigenvalue is e1 = e[1] and the first 
# eigenvector is v̄1 = U̲[:,1]

# We can check that the secular equation |A̲ - eI̲| = 0
# is satisfied (for each eigenvalue)
for i in 1:4
   isapprox( det(A - e[i]*I), 0.00, atol=1e-4) || error("Secular equation $i not satisfied")
end

# The corresponding eigenvectors satisfy (A̲ - eI̲)v̄ = 0̄
for i in 1:4
   isapprox( (A - e[i]*I)*U[:,i], fill(0.00,4), atol=1e-4) || error("a) Eigenvalue equation $i not satisfied")
end
# or equivalently A̲v̄ = ev̄
for i in 1:4
   isapprox( A*U[:,i], e[i]*U[:,i], atol=1e-4) || error("b) Eigenvalue equation $i not satisfied")
end

# - If A̲ is symmetric, thus U̲ is orthonormal and U̲' = U̲^-1
# - If A̲ is hermitic, thus U̲ is unitary and U̲' = U̲^-1
#
# In both cases it can be described as a similarity transformation
#  A̲ = U̲ D̲ U̲^-1 = U̲ D̲ U̲'
# where D is a diagonal matrix with the eigenvalues as diagonal elements.
# Therefore, U is the matrix that diagonalizes A.
# A̲ and D̲ are said to be similar matrices.
#
A ≈ U * Diagonal(e) * U' || error("Not a similarity transformation")

# A similarity transformation is also a conformal mapping. A matrix W̲
# encoding a crystallographic symmetry operation is an isometry if
# G̲ ≈ W̲ G̲ W̲'
# where G̲ is the metric tensor.





