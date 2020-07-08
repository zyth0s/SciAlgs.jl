
# Helium ground state
# The most naive implementation
# - Does not take advantage of ERIs symmetry
# - Fixed number of SCF steps without checking convergence
# - Unoptimized
# - Simple basis set
# Computational Physics by Jos Thijssen, §4.3.2

import LinearAlgebra: eigen

# Basis set
α1 =  0.297104
α2 =  1.236745
α3 =  5.749982
α4 = 38.216677
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
  -4π/(αp+αq)
end

function eri_over_primitives(αp,αr,αq,αs)
  2π^(5/2)/((αp+αq)*(αr+αs)*sqrt(αp+αq+αr+αs))
end


# Compute one-electron matrices

h=zeros(nbasis,nbasis)
S=zeros(nbasis,nbasis)
Q=zeros(nbasis,nbasis,nbasis,nbasis)

for (p,αp) in enumerate(exps), 
    (q,αq) in enumerate(exps)
  S[p,q] = overlap_over_primitives(αp,αq)
  h[p,q] = kinetic_over_primitives(αp,αq) + nuclear_attraction_over_primitives(αp,αq)
end

# Compute two-electron interaction tensor
energies,C = eigen(h,S)

for (p,αp) in enumerate(exps), 
    (r,αr) in enumerate(exps),
    (q,αq) in enumerate(exps),
    (s,αs) in enumerate(exps)
  Q[p,r,q,s] = eri_over_primitives(αp,αr,αq,αs)
end

# SCF
function SCF(S,h,Q,C)
  nbasis = size(h)[1]
  for _ in 1:20 # Most naive approach
    F = zeros(nbasis,nbasis)
    for p in 1:nbasis, 
        q in 1:nbasis
        F[p,q] = h[p,q]
        for r in 1:nbasis,
            s in 1:nbasis
          F[p,q] += Q[p,r,q,s]*C[r,1]*C[s,1]
        end
    end
    energies,C = eigen(F,S)
  end
  energies,C,F
end

energies,C,F = SCF(S,h,Q,C)

# Ground state energy
function ground_energy(h,Q,C)
  nbasis = size(h)[1]
  E_G = 0.0
  for p in 1:nbasis, 
      q in 1:nbasis
      E_G += 2C[p,1]*C[q,1]*h[p,q]
  end
  for p in 1:nbasis, 
      q in 1:nbasis,
      r in 1:nbasis,
      s in 1:nbasis
    E_G += Q[p,r,q,s]*C[p,1]*C[q,1]*C[r,1]*C[s,1]
  end
  E_G
end

Eerror = ground_energy(h,Q,C) - -2.8551716
# Meeh ... accuracy
isapprox(ground_energy(h,Q,C),-2.8551716,atol=2e-5) || error("Wrong energy") # Hartree




