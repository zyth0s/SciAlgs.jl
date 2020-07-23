
# Calculation of the electronic structure of a molecule
# with the help of PySCF to calculate AO integrals and hold molecular data
# * Restricted Hartree-Fock
# * Møller-Plesset order 2
# * Coupled Cluster Singles and Doubles
# Calculates the energy of any molecule.
# Instructions at https://github.com/CrawfordGroup/ProgrammingProjects
# may be of interest.

using PyCall: pyimport
using Formatting: printfmt
using LinearAlgebra: Diagonal, Hermitian, eigen, norm, tr, diag, dot

pyscf = pyimport("pyscf")
mp = pyimport("pyscf.mp")   # Had to import mp alone ??!
cc = pyimport("pyscf.cc") # Had to import mp alone ??!
np = pyimport("numpy") # alternative: Einsum.jl

function pyscf_atom_from_xyz(fpath::String)
   join(split(read(open(fpath),String),"\n")[3:end],"\n")
end

function index(i,j)
   m,M = minmax(i-1,j-1)
   M*(M+1)÷2 + m + 1
end

# Compound indices ijkl
function get_4idx(i,j,k,l)
   ij = index(i,j)
   kl = index(k,l)
   index(ij,kl)
end

# Löwdin symmetric Orthogonalization (with A = X^(-1/2))
function lowdinOrtho(A,B)
   A' * B * A
end

function buildFock(hcore,D,nao,eri)
   F = copy(hcore)
   for i in 1:nao, j in 1:nao,
       k in 1:nao, l in 1:nao
      ijkl = get_4idx(i,j,k,l)
      ikjl = get_4idx(i,k,j,l)
      F[i,j] += D[k,l]*(2eri[ijkl]-eri[ikjl])
   end
   F
end

#########################################################################
# Project #3: The Hartree-Fock self-consistent field (SCF) procedure.
#########################################################################
function scf_rhf(mol)
   nelec = mol.nelectron
   nocc  = mol.nelectron ÷ 2 
   nao   = mol.nao_nr() |> n -> convert(Int,n) # why not converted ??!
   enuc  = mol.energy_nuc()
   s     = mol.intor("cint1e_ovlp_sph")
   hcore = mol.intor("cint1e_nuc_sph") + mol.intor("cint1e_kin_sph")
   eri   = mol.intor("cint2e_sph", aosym="s8")

   ########################################################
   #4: BUILD THE ORTHOGONALIZATION MATRIX
   ########################################################
   # NOTE: The AO basis is non-orthogonal

   # Symmetric orthogonalization matrix
   s_minushalf = s^(-0.5) # Matrix power (via JordanNF)

   ########################################################
   #5: BUILD THE INITIAL (GUESS) DENSITY
   ########################################################

   P = zeros(nao,nao) # D = 2P, P often called density matrix

   ########################################################
   #6: COMPUTE THE INITIAL SCF ENERGY
   ########################################################

   Eelec = 0.0

   #*******************************************************
   #*******************************************************
   #
   #             SCF ITERATIVE PROCESS
   #
   #*******************************************************
   #*******************************************************

   iteri     = 0    # extant iteration
   itermax   = 100  # max number of iterations
   δE        = 1e-8 # energy threshold for convergence
   δP        = 1e-8 # P-matrix threshold for convergence
   ΔE        = 1.0; @assert ΔE > δE # initial E difference
   ΔP        = 1.0; @assert ΔP > δP # initial P difference
   Eelec     = 0.0
   oldEelec  = 0.0
   Etot      = 0.0  # total energy
   e         = zeros(nao)
   C         = zeros(nao,nao)
   DIISstack = [] # for DIIS
   Fstack    = [] # for DIIS
   println(" Iter        E(elec)        E(tot)        ΔE(elec)         ΔP")
   println("-----   --------------  -------------  ------------  -------------")

   while abs(ΔE) > δE && iteri < itermax && ΔP > δP

      iteri += 1

      ########################################################
      #7: COMPUTE THE NEW FOCK MATRIX
      ########################################################

      F = buildFock(hcore,P,nao,eri)

      #########################################################################
      # Project #8: DIIS extrapolation for the SCF procedure.
      # - P. Pulay, Chem. Phys. Lett. 73, 393-398 (1980).
      # - P. Pulay, J. Comp. Chem. 3, 556-560 (1982).
      # - T. P. Hamilton and P. Pulay, J. Chem. Phys. 84, 5728-5734 (1986).
      # - C. David Sherrill. "Some comments on accellerating convergence of iterative
      #      sequences using direct inversion of the iterative subspace (DIIS)".
      #      Available at: vergil.chemistry.gatech.edu/notes/diis/diis.pdf. (1998)
      #########################################################################
      # orbital gradient in AO basis: F Di Si - S Di Fi
      # better choice:     (S^-0.5)' (F Di Si - S Di Fi) S^-0.5
      DIIS_residual = s * P * F
      DIIS_residual = DIIS_residual' - DIIS_residual
      DIIS_residual = s_minushalf * DIIS_residual * s_minushalf
      if iteri > 1
         push!(Fstack,F)
         push!(DIISstack, DIIS_residual)
      end
      ΔDIIS = norm(DIIS_residual)

      if iteri > 2
         if iteri > 15 # Limit the storage of DIIS arrays
            popfirst!(DIISstack)
            popfirst!(Fstack)
         end
         dim_B = length(Fstack) + 1 # 1 dim for the Lagrange multiplier
         B = zeros(dim_B,dim_B)
         B[end,:]   .= -1
         B[:,  end] .= -1
         B[end,end]  =  0
         for i in eachindex(Fstack), j in eachindex(Fstack)
            B[i,j] = dot(DIISstack[i], DIISstack[j]) # Gramian matrix
         end
         # Solve Lagrange equation of Pulay
         Pulay_rhs = zeros(dim_B)
         Pulay_rhs[end] = -1
         coef_Pulay = B \ Pulay_rhs
         F = zeros(size(F))
         for (i,c) in enumerate(coef_Pulay[1:end-1]) # skip Lagrange mult.
            F += c * Fstack[i]
         end
      end

      ########################################################
      #8: BUILD THE NEW DENSITY MATRIX
      ########################################################

      # Transformed Fock matrix
      Fp = lowdinOrtho(s_minushalf,F)
      # MO coeffs as linear combination of orthonormal MO basis functions
      e, Cp = eigen(Hermitian(Fp))
      # New MO Coefficients as linear combination of AO basis functions
      C = s_minushalf * Cp
      # Renew P-matrix.
      P, oldP = zeros(nao,nao), copy(P)

      mocc = C[:,1:nocc]
      P  = mocc * mocc'
      ΔP = norm(P - oldP) # Frobenius norm

      ########################################################
      #9: COMPUTE THE NEW SCF ENERGY
      ########################################################

      Eelec, oldEelec = tr(P * (hcore + F)), Eelec # eq. (238) Janos

      ΔE   = Eelec - oldEelec
      Etot = Eelec + enuc

      printfmt(" {1:3d}  {2:14.8f}  {3:14.8f}  {4:14.8f} {5:14.8f}\n",
                 iteri,   Eelec,       Etot,      ΔE   ,    ΔP)
      ########################################################
      #10: TEST FOR CONVERGENCE: WHILE LOOP (go up)
      ########################################################
   end
   @assert iteri < itermax "NOT CONVERGED!!"

   dipole = dipole_moment(2P,mol)
   printfmt("The dipole moment (Debye): {1:8.4f} {2:8.4f} {3:8.4f}\n",dipole...)

   # In RHF, D = 2P, an occupied spatial orbital means 2 elec, α and β
   # PS is a "correct" representation of the 1RDM in AO basis
   Mulliken = 1; Lowdin = 0.5
   println("Mulliken population analysis")
   population_analysis(2P,s,Mulliken,mol)
   println("Löwdin population analysis")
   population_analysis(2P,s,Lowdin,mol)

   Etot,e,C
end

#########################################################################
# Additional Concepts
# also see: "Hartree-Fock methods of Quantum Chemistry" by Janos Angyan
#########################################################################
# One-electron properties
function dipole_moment(D, mol)
   # eq. (365) Janos
   # It does the same as PySCF does.
   #ao_dip = @pywith mol.with_common_orig([0,0,0]) begin
   #   mol.intor_symmetric("int1e_r", comp=3)
   #end
   common_orig_bak = mol._env[1:3]
   mol = mol.set_common_orig([0,0,0])
   ao_dip = mol.intor_symmetric("int1e_r", comp=3)
   mol = mol.set_common_orig(common_orig_bak)

   el_dip   = real.(np.einsum("xij,ji->x", ao_dip, D))
   charges  = mol.atom_charges()
   coords   = mol.atom_coords()
   nucl_dip = np.einsum("i,ix->x", charges, coords)
   mol_dip  = nucl_dip - el_dip
   mol_dip *= pyscf.data.nist.AU2DEBYE
end

# Population Analysis/Atomic Charges
@doc raw"""
   `population_analysis(D,S,α,mol)`

Generalized Mulliken, Lowdin, ... population analysis
with `D` the AO density matrix, `S` the overlap matrix,
`α` defines the analysis type, and mol is a PySCF molecule.

Qa = Za - ∑μ∈a [S^α P S^(1-α)]μμ  with 0 < α < 1.

In particular:
* α = 1 => Mulliken
* α = ½ => Löwdin
"""
function population_analysis(D,S,α,mol)
   s_α, s_mα = S^α, S^(1-α) # Matrix power (via JordanNF)
   @assert isapprox(mol.nelectron, tr(s_α*D*s_mα), atol=1e-8) # eq. (363) Janos
   PopAO = diag(s_α * D * s_mα)
   Q = zeros(mol.natm)
   for (μ,record) in enumerate(mol.ao_labels(fmt=nothing))
      a = record[1] + 1 |> n -> convert(Int,n)
      Q[a] -= PopAO[μ] # eqs. (362,364) Janos
   end
   Q += mol.atom_charges()
   for a in eachindex(Q)
      println(" Atom $a has charge $(Q[a])")
   end
   @assert isapprox(sum(Q), mol.charge, atol=1e-8)
end

#########################################################################
# Project #4: The Second-Order Moller-Plesset
#             Perturbation Theory (MP₂) Energy
#########################################################################

########################################################
#3: TRANSFORM THE TWO-ELECTRON INTEGRALS TO THE MO BASIS
########################################################

function ao2mo_noddy(nao,C,eri,mol)
   eri_mo = zeros(length(eri))
   @inbounds for q in 1:nao, p in q:nao
      @inbounds for r in 1:p
         lim = ifelse(p == r, q, r)
         @inbounds for s in 1:lim
            pqrs = get_4idx(p,q,r,s)
            @inbounds for i in 1:nao, j in 1:nao, k in 1:nao, l in 1:nao
               ijkl = get_4idx(i,j,k,l)
               eri_mo[pqrs] += C[i,p]*C[j,q]*eri[ijkl]*C[k,r]*C[l,s]
            end
         end
      end
   end
   eri_mo
end

function ao2mo_smart(nao,C,eri,mol)
   """
   (ij|kl) => (ij|ks) => (ij|rs) => (iq|rs) => (pq|rs)
   """
   M      = nao*(nao+1)÷2
   X = zeros(nao, nao)
   tmp = zeros(M,M)
   @inbounds for j in 1:nao, i in j:nao
      ij = index(i,j)
      @inbounds for l in 1:nao, k in l:nao
         ijkl = get_4idx(i,j,k,l)
         X[k,l] = X[l,k] = eri[ijkl]
      end

      Y = zeros(nao, nao)
      Y = C' * X
      X = zeros(nao, nao)
      X = Y * C
      @inbounds for l in 1:nao, k in l:nao
         kl = index(k,l)
         tmp[kl,ij] = X[k,l]
      end
   end

   eri_mo = zeros(M*(M+1)÷2)

   @inbounds for l in 1:nao, k in l:nao
      kl = index(k,l)
      X = zeros(nao,nao)
      Y = zeros(nao,nao)
      @inbounds for j in 1:nao, i in j:nao
         ij = index(i,j)
         X[i,j] = X[j,i] = tmp[kl,ij]
      end
      Y = zeros(nao, nao)
      Y = C' *  X
      X = zeros(nao, nao)
      X = Y * C
      @inbounds for j in 1:nao, i in j:nao
         klij = get_4idx(k,l,i,j)
         eri_mo[klij] = X[i,j]
      end
   end
   eri_mo
end
#eri_mo   = pyscf.ao2mo(nao,C,eri,mol)

########################################################
#4: COMPUTE THE MP₂ ENERGY
########################################################
function short_mp2(e,C,eri,mol)
   nocc  = mol.nelectron ÷ 2
   nao   = mol.nao_nr() |> n -> convert(Int,n) # why not converted ??!
   # Convert AO -> MO basis
   eri_mo = zeros(nao,nocc,nao,nocc)
   emp2 = 0.0
   occupied =      1:nocc
   virtual  = nocc+1:nao
   eri = pyscf.ao2mo.restore(1,eri,nao) # reshape
   eri_mo = np.einsum("pa,qi,rb,sj,pqrs->aibj",C[:,1:nocc],C,C[:,1:nocc],C,eri)
   @inbounds for i in occupied, a in virtual,
                 j in occupied, b in virtual
       eiajb, eibja = eri_mo[i,a,j,b],  eri_mo[i,b,j,a]
       emp2 += eiajb * (2.0eiajb-eibja) / (e[i]+e[j]-e[a]-e[b])
    end
   emp2
end

function mymp2(e,C,eri,mol; alg::Symbol=:ao2mo_smart)
   # alg ∈ [ao2mo_noddy, ao2mo_smart]
   ao2mo    = getfield(Main, alg) # get function
   nocc     = mol.nelectron ÷ 2
   nao      = mol.nao_nr()       |> n -> convert(Int,n) # why not converted ??!
   emp2     = 0.0
   occupied =      1:nocc
   virtual  = nocc+1:nao
   eri_mo   = ao2mo(nao,C,eri,mol)
   @inbounds for i in occupied, a in virtual,
                 j in occupied, b in virtual
      iajb = get_4idx(i,a,j,b)
      ibja = get_4idx(i,b,j,a)
      emp2 += eri_mo[iajb] * (2.0eri_mo[iajb]-eri_mo[ibja]) / (e[i]+e[j]-e[a]-e[b])
   end

   # In the basis of spin-molecular orbitals:
   #emp2 += 0.25eri_mo[ijab]*eri_mo[ijab]/(e[i]+e[j]-e[a]-e[b])
   emp2
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Project #5: The Coupled Cluster Singles and Doubles (CCSD) Energy
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    J.F. Stanton, J. Gauss, J.D. Watts,
#    and R.J. Bartlett, J. Chem. Phys.
#    volume 94, pp. 4334-4345 (1991).

# We follow the convention that
# i, j, k,... represent occupied orbitals, with
# a, b, c,... unoccupied.
# p, q, r,... are generic indices which may represent
#             either occupied or unoccupied
# orbitals.

###################################################
#1: Preparing the Spin-Orbital Basis Integrals
###################################################
# Translate integrals from spatial mo basis to spin-orbital basis
function orb_to_spinorb(e,eri_mo)
   # 1  2  3  4  5  6  7  8
   # α, β, α, β, α, β, α, β, ...
   nao  = length(e)
   nsmo = 2nao
   eri_smo = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for p in 1:nsmo, q in 1:nsmo, r in 1:nsmo, s in 1:nsmo
      prqs = get_4idx((p+1)÷2,(r+1)÷2,(q+1)÷2,(s+1)÷2)
      psqr = get_4idx((p+1)÷2,(s+1)÷2,(q+1)÷2,(r+1)÷2)
      value1 = eri_mo[prqs] * (mod(p,2) == mod(r,2)) * (mod(q,2) == mod(s,2))
      value2 = eri_mo[psqr] * (mod(p,2) == mod(s,2)) * (mod(q,2) == mod(r,2))
      eri_smo[p,q,r,s] = value1 - value2
   end
   # The fock matrix is diagonal in the spin-MO basis
   fs = zeros(nsmo)
   for i in 1:nsmo
      fs[i] = e[(i+1)÷2]
   end
   fs = Diagonal(fs)
   eri_smo, fs
end

###################################################
#2: Build the Initial-Guess Cluster Amplitudes
###################################################

function initial_amplitudes(occupied,virtual,fs,eri_smo)
   nsmo = length(occupied) + length(virtual)
   t1 = zeros(nsmo,nsmo)
   t2 = zeros(nsmo,nsmo,nsmo,nsmo)

   # MP₂ guess for T₂ amplitudes
   # t[a,b,i,j] = ⟨ij||ab⟩/ (ϵi + ϵj - ϵa - ϵb)
   # E_mp2      = 1/4 * ∑ijab ⟨ij||ab⟩ t[a,b,i,j]
   emp2 = 0.0
   @inbounds for i in occupied, a in virtual,
                 j in occupied, b in virtual
      t2[a,b,i,j] += eri_smo[i,j,a,b] / (fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b])
      emp2 += 0.5*0.5*t2[a,b,i,j]*eri_smo[i,j,a,b]
   end
   t1, t2, emp2
end

###################################################
#3: Calculate the CC Intermediates
###################################################

# Effective doubles
function tau_(t1,t2,a,b,i,j)
   # Stanton1991 (9)
   t2[a,b,i,j] + 0.5(t1[a,i]*t1[b,j] - t1[b,i]*t1[a,j])
end
function tau(t1,t2,a,b,i,j)
   # Stanton1991 (10)
   t2[a,b,i,j] + t1[a,i]*t1[b,j] - t1[b,i]*t1[a,j]
end

function updateF_W(occupied,virtual, # occ and virt spin-MO spaces
                     t1,t2,            # CC amplitudes
                     fs,eri_smo)       # integrals in spin-MO basis
   nsmo = length(occupied) + length(virtual)
   # Stanton1991 (3)
   Fae = zeros(nsmo,nsmo)
   @inbounds for a in virtual, e in virtual
      Fae[a,e] = (1 - (a == e))*fs[a,e]
      @inbounds for m in occupied
         Fae[a,e] += -0.5*fs[m,e]*t1[a,m]
         @inbounds for f in virtual
            Fae[a,e] += t1[f,m]*eri_smo[m,a,f,e]
            @inbounds for n in occupied
               Fae[a,e] += -0.5*tau_(t1,t2,a,f,m,n)*eri_smo[m,n,e,f]
            end
         end
      end
   end
   # Stanton1991 (4)
   Fmi = zeros(nsmo,nsmo)
   @inbounds for m in occupied, i in occupied
      Fmi[m,i] = (1 - (m == i))*fs[m,i]
      @inbounds for e in virtual
         Fmi[m,i] += 0.5t1[e,i]*fs[m,e]
         @inbounds for n in occupied
            Fmi[m,i] += t1[e,n]*eri_smo[m,n,i,e]
            @inbounds for f in virtual
               Fmi[m,i] += 0.5tau_(t1,t2,e,f,i,n)*eri_smo[m,n,e,f]
            end
         end
      end
   end
   # Stanton1991 (5)
   Fme = zeros(nsmo,nsmo)
   @inbounds for e in virtual, m in occupied
      Fme[m,e] = fs[m,e]
      @inbounds for f in virtual, n in occupied
         Fme[m,e] += t1[f,n]*eri_smo[m,n,e,f]
      end
   end
   # Stanton1991 (6)
   Wmnij = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for m in occupied, n in occupied, i in occupied, j in occupied
      Wmnij[m,n,i,j] = eri_smo[m,n,i,j]
      @inbounds for e in virtual # P_(ij)
         Wmnij[m,n,i,j] += t1[e,j]*eri_smo[m,n,i,e] -
                           t1[e,i]*eri_smo[m,n,j,e]
         @inbounds for f in virtual
            Wmnij[m,n,i,j] += 0.25tau(t1,t2,e,f,i,j)*eri_smo[m,n,e,f]
         end
      end
   end
   # Stanton1991 (7)
   Wabef = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for a in virtual, b in virtual, e in virtual, f in virtual
      Wabef[a,b,e,f] = eri_smo[a,b,e,f]
      @inbounds for m in occupied # P_(ab)
         Wabef[a,b,e,f] += -t1[b,m]*eri_smo[a,m,e,f] +
                            t1[a,m]*eri_smo[b,m,e,f]
         @inbounds for n in occupied
            Wabef[a,b,e,f] += 0.25tau(t1,t2,a,b,m,n)*eri_smo[m,n,e,f]
         end
      end
   end
   # Stanton1991 (8)
   Wmbej = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for b in virtual, m in occupied, e in virtual, j in occupied
      Wmbej[m,b,e,j] = eri_smo[m,b,e,j]
      @inbounds for f in virtual
         Wmbej[m,b,e,j] += t1[f,j]*eri_smo[m,b,e,f]
      end
      @inbounds for n in occupied
         Wmbej[m,b,e,j] += -t1[b,n]*eri_smo[m,n,e,j]
         @inbounds for f in virtual
            Wmbej[m,b,e,j] += -(0.5*t2[f,b,j,n] + t1[f,j]*t1[b,n])*eri_smo[m,n,e,f]
         end
      end
   end
   Fae, Fmi, Fme, Wmnij, Wabef, Wmbej
end

###############################################
#4: Compute the Updated Cluster Amplitudes
###############################################

# Stanton1991 (1)
function updateT1(occupied,virtual, # occ and virt spin-MO spaces
                  t1,t2,            # CC amplitudes
                  fs,eri_smo,       # integrals in spin-MO basis
                  Fae,Fmi,Fme)      # one-particle CC intermediates
   nsmo = length(occupied) + length(virtual)
   # Denominator energies
   Dai = zeros(nsmo,nsmo)
   @inbounds for a in virtual,i in occupied
      # Stanton1991 (12)
      Dai[a,i] = fs[i,i] - fs[a,a]
   end
   _t1 = zeros(nsmo,nsmo)
   @inbounds for a in virtual, i in occupied
      _t1[a,i] = fs[i,a] # 1st RHS term (leading term in the expansion)
      @inbounds for e in virtual
         _t1[a,i] += t1[e,i]*Fae[a,e] # 2nd RHS term
      end
      @inbounds for m in occupied
         _t1[a,i] += -t1[a,m]*Fmi[m,i] # 3rd RHS term
         @inbounds for e in virtual
            _t1[a,i] += t2[a,e,i,m]*Fme[m,e] # 4th RHS term
            @inbounds for f in virtual
               _t1[a,i] += -0.5t2[e,f,i,m]*eri_smo[m,a,e,f] # 6th RHS term
            end
            @inbounds for n in occupied
               _t1[a,i] += -0.5t2[a,e,m,n]*eri_smo[n,m,e,i] # 7th RHS term
            end
         end
      end
      @inbounds for f in virtual,n in occupied
         _t1[a,i] += -t1[f,n]*eri_smo[n,a,i,f] # 5th RHS term
      end
      _t1[a,i] /= Dai[a,i] # LHS
   end
   _t1
end

# Stanton1991 (2)
function updateT2(occupied,virtual,  # occ and virt spin-MO spaces
                  t1,t2,             # CC amplitudes
                  fs,eri_smo,        # integrals in spin-MO basis
                  Fae,Fmi,Fme,       # one-particle CC intermediates
                  Wmnij,Wabef,Wmbej) # two-particle CC intermediates
   nsmo = length(occupied) + length(virtual)
   # Denominator energies
   Dabij = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for a in virtual, i in occupied, b in virtual, j in occupied
      # Stanton1991 (13)
      Dabij[a,b,i,j] = fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b]
   end
   _t2 = zeros(nsmo,nsmo,nsmo,nsmo)
   @inbounds for a in virtual, i in occupied, b in virtual, j in occupied
      _t2[a,b,i,j] = eri_smo[i,j,a,b] # 1st RHS (leading term in the expansion)
      @inbounds for e in virtual # 2nd RHS term
         taeij, tbeij = t2[a,e,i,j], t2[b,e,i,j]
         _t2[a,b,i,j] += taeij*Fae[b,e] -
                         tbeij*Fae[a,e]   # asym. permutes b <-> a
         @inbounds for m in occupied # 3rd RHS term
            _t2[a,b,i,j] += -0.5taeij*t1[b,m]*Fme[m,e] +
                             0.5tbeij*t1[a,m]*Fme[m,e]   # asym. permutes b <-> a
         end
      end
      @inbounds for m in occupied # 4th RHS term
         tabim, tabjm = t2[a,b,i,m], t2[a,b,j,m]
         _t2[a,b,i,j] += -tabim*Fmi[m,j] +
                          tabjm*Fmi[m,i]  # asym. permutes i <-> j
         @inbounds for e in virtual # 5th RHS term
            _t2[a,b,i,j] += -0.5tabim*t1[e,j]*Fme[m,e] +
                             0.5tabjm*t1[e,i]*Fme[m,e]  # asym. permutes i <-> -j
         end
      end
      @inbounds for e in virtual # 10th RHS term
         _t2[a,b,i,j] += t1[e,i]*eri_smo[a,b,e,j] -
                         t1[e,j]*eri_smo[a,b,e,i]  # asym. permutes i <-> j
         @inbounds for f in virtual # 7th RHS term
            _t2[a,b,i,j] += 0.5tau(t1,t2,e,f,i,j)*Wabef[a,b,e,f]
         end
      end
      @inbounds for m in occupied # 11th RHS term
         _t2[a,b,i,j] += -t1[a,m]*eri_smo[m,b,i,j] +
                          t1[b,m]*eri_smo[m,a,i,j]   # asym. permutes a <-> b
         @inbounds for e in virtual # 8-9th RHS terms
            taeim, taejm = t2[a,e,i,m], t2[a,e,j,m]
            tbeim, tbejm = t2[b,e,i,m], t2[b,e,j,m]
            _t2[a,b,i,j] +=  taeim*Wmbej[m,b,e,j] - t1[e,i]*t1[a,m]*eri_smo[m,b,e,j] +
                            -taejm*Wmbej[m,b,e,i] + t1[e,j]*t1[a,m]*eri_smo[m,b,e,i] + # i <-> j
                            -tbeim*Wmbej[m,a,e,j] - t1[e,i]*t1[b,m]*eri_smo[m,a,e,j] + #          a <-> b
                             tbejm*Wmbej[m,a,e,i] - t1[e,j]*t1[b,m]*eri_smo[m,a,e,i]   # i <-> j, a <-> b
         end
         @inbounds for n in occupied # 6th RHS term
            _t2[a,b,i,j] += 0.5tau(t1,t2,a,b,m,n)*Wmnij[m,n,i,j]
         end
      end
      _t2[a,b,i,j] /= Dabij[a,b,i,j] # LHS
   end
   _t2
end

###########################################
#5: Check for Convergence and Iterate
###########################################

function extant_E_CCSD(occupied, virtual, fs, t1, t2, eri_smo)
   # Calculate the extant CC correlation energy
   # E_CC = ∑ia fia t[a,i] + 1/4 ∑ijab ⟨ij||ab⟩t[a,b,i,j] + 1/2 ∑ijab ⟨ij||ab⟩ t[a,i] t[b,j]
   E_CCSD = 0.0
   @inbounds for a in virtual, i in occupied
      E_CCSD += fs[i,a] * t1[a,i]
      @inbounds for b in virtual, j in occupied
         E_CCSD += 0.25eri_smo[i,j,a,b] * t2[a,b,i,j] +
                    0.5eri_smo[i,j,a,b] * t1[a,i] * t1[b,j]
      end
   end
   E_CCSD
end

function ccsd(e, C, eri)

   nao      = size(C)[1]
   nsmo     = 2nao
   nelec    = mol.nelectron
   occupied =       1:nelec
   virtual  = nelec+1:nsmo
   eri_mo = ao2mo_smart(nao,C,eri,mol)

   eri_smo, fs = orb_to_spinorb(e,eri_mo)

   t1,t2, E_MP2 = initial_amplitudes(occupied,virtual,fs,eri_smo)

   E_CCSD  = 0.0
   ΔE_CCSD = 1.0
   iteri   = 0
   itermax = 60
   println(" Iter        E(CCSD)        ΔE(CCSD)")
   println("-----   --------------   --------------")
   while ΔE_CCSD > 1e-9 && iteri < itermax
      iteri += 1
      Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateF_W(occupied,virtual,t1,t2,fs,eri_smo)
      t1 = updateT1(occupied,virtual,t1,t2,fs,eri_smo,Fae,Fmi,Fme)
      t2 = updateT2(occupied,virtual,t1,t2,fs,eri_smo,Fae,Fmi,Fme,Wmnij,Wabef,Wmbej)
      E_CCSD, oldE_CCSD = extant_E_CCSD(occupied,virtual,fs,t1,t2,eri_smo), E_CCSD
      ΔE_CCSD = abs(E_CCSD - oldE_CCSD)
      printfmt("{1:5d}   {2:13.10f}   {3:13.10f}\n",iteri,E_CCSD,ΔE_CCSD)
   end
   @assert iteri < itermax "NOT CONVERGED!!"
   E_CCSD, E_MP2
end


#########################################################################
# Test HF, MP₂, and CCSD
mol = pyscf.gto.Mole()
mol.build(verbose = 0,
          basis = "sto-3g",
          #atom = pyscf_atom_from_xyz("../data/h2.xyz"),
          #atom = pyscf_atom_from_xyz("../data/acetaldehyde.xyz"),
          atom = """8 0  0.    0
                    1 0 -0.757 0.587
                    1 0  0.757 0.587""",
         )

mf = pyscf.scf.RHF(mol)

e_rhf_me, e, C = scf_rhf(mol)
e_rhf_pyscf    = mf.kernel()

mf.analyze(verbose=3,ncol=10, digits=9)

println( "           E(HF)     ")
println( "        -------------")
printfmt(" ME:    {1:13.8f}  \n", e_rhf_me)
printfmt(" PySCF: {1:13.8f}  \n", e_rhf_pyscf)

eri   = mol.intor("cint2e_sph", aosym="s8")
println()
e_mp2_me1   = @time mymp2(e,C,eri,mol,alg=:ao2mo_noddy) # Naive  ao2mo with my MP₂
e_mp2_me2   = @time mymp2(e,C,eri,mol,alg=:ao2mo_smart) # Smart  ao2mo with my MP₂
e_mp2_me3   = @time short_mp2(e,C,eri,mol)              # einsum ao2mo with my MP₂
e_mp2_pyscf = @time mp.MP2(mf).kernel()[1]              # pure PySCF MP₂ solution
println()

println( "             ΔE(MP₂)        E(MP₂)   ")
println( "          -----------   -------------")
printfmt(" ME1:   {1:13.8f}   {2:13.8f}  \n", e_mp2_me1,  e_rhf_me   +e_mp2_me1)
printfmt(" ME2:   {1:13.8f}   {2:13.8f}  \n", e_mp2_me2,  e_rhf_me   +e_mp2_me1)
printfmt(" ME3:   {1:13.8f}   {2:13.8f}  \n", e_mp2_me3,  e_rhf_me   +e_mp2_me1)
printfmt(" PySCF: {1:13.8f}   {2:13.8f}  \n", e_mp2_pyscf,e_rhf_pyscf+e_mp2_pyscf)
println()

@time e_ccsd_me, emp2 = ccsd(e,C,eri)
@time e_ccsd_pyscf    = cc.ccsd.CC(mf).kernel()[1]
emp2 ≈ e_mp2_me3 || @warn "E_MP2 from CC is not accurate"
println()

println("            E(HF)            ΔE(MP2)           E(MP2)           ΔE(CCSD)  ")
println("        --------------   --------------   ---------------   ---------------")
printfmt(" ME:    {1:13.10f}    {2:13.10f}    {3:13.10f}    {4:13.10f}\n",
                 e_rhf_me,      e_mp2_me3,     e_rhf_me+e_mp2_me3,    e_ccsd_me)
printfmt(" PySCF: {1:13.10f}    {2:13.10f}    {3:13.10f}    {4:13.10f}\n",
                  e_rhf_pyscf,  e_mp2_pyscf,e_rhf_pyscf+e_mp2_pyscf,e_ccsd_pyscf)
