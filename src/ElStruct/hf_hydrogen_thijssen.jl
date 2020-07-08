
# Hydrogen ground state
# The most naive implementation
# - Unoptimized
# - Simple basis set
# Computational Physics by Jos Thijssen, §3.2.2

import LinearAlgebra: eigen

# Basis set
α1 = 13.00773
α2 =  1.962079
α3 =  0.444529
α4 =  0.1219492
exps = [α1,α2,α3,α4]
nbasis = length(exps)

function basis_function(α,r)
  exp(-α*r^2)
end

# Integrals over primitives
function overlap_over_primitives(αp,αq)
  (π/(αp+αq))^(3/2)
end

function kinetic_over_primitives(αp,αq)
  3αp*αq*π^(3/2)/(αp+αq)^(5/2)
end

function nuclear_attraction_over_primitives(αp,αq)
  -2π/(αp+αq)
end

# Compute one-electron matrices

h=zeros(nbasis,nbasis)
S=zeros(nbasis,nbasis)

for (p,αp) in enumerate(exps), 
    (q,αq) in enumerate(exps)
  S[p,q] = overlap_over_primitives(αp,αq)
  h[p,q] = kinetic_over_primitives(αp,αq) + nuclear_attraction_over_primitives(αp,αq)
end

# Diagonalize the Fock matrix, F = h
energies,C = eigen(h,S)

# Test
Eerror = energies[1] - -0.499278
isapprox(energies[1],-0.499278,atol=1e-6) || error("Wrong energy") # Hartree




