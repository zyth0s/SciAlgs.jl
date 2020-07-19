
# Restricted Hartree-Fock program that uses PySCF to calculate the integrals in the AO basis.
# Calculates the energy of any molecule, however, it has no SCF speedup.
# Therefore, it may converge poorly.
# Instructions at https://github.com/zyth0s/ProgrammingProjects/tree/master/Project%2303
# may be of interest.

using PyCall: pyimport
using Formatting: printfmt
using LinearAlgebra: Diagonal, Hermitian, eigen, norm

pyscf = pyimport("pyscf")
mp = pyimport("pyscf.mp") # Had to import mp alone ??!
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

# Lowdin symmetric Orthogonalization (with A = X^(-1/2)) 
function lowdinOrtho(A,B)
   A' * B * A
end

function buildFock(hcore,D,nao,eri)
   F = deepcopy(hcore)
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

   # Diagonalize the overlap matrix
   s_diag, s_eigvec = eigen(s)
   s_diag_minushalf = Diagonal(s_diag.^(-0.5))
   # Build the symmetric orthogonalization matrix
   s_minushalf = s_eigvec * s_diag_minushalf * s_eigvec'

   ########################################################
   #5: BUILD THE INITIAL (GUESS) DENSITY
   ########################################################

   D = zeros(nao,nao)

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

   iteri     = 0    # current iteration
   itermax   = 100  # max number of iterations
   δE        = 1e-8 # energy threshold for convergence
   δD        = 1e-8 # density matrix (D) threshold for convergence
   ΔE        = 1.0; @assert ΔE > δE # initial E difference
   ΔD        = 1.0; @assert ΔD > δD # initial D difference
   Etot      = 0.0  # total energy
   e         = zeros(nao)
   C         = zeros(nao,nao)
   println(" Iter        E(elec)        E(tot)        Delta(E)        RMS(D)")
   println("-----   --------------  -------------  ------------  -------------")

   while abs(ΔE) > δE && iteri < itermax && ΔD > δD

      iteri += 1

      ########################################################
      #7: COMPUTE THE NEW FOCK MATRIX
      ########################################################

      F = buildFock(hcore,D,nao,eri)

      ########################################################
      #8: BUILD THE NEW DENSITY MATRIX
      ########################################################

      # Transformed Fock matrix
      Fp = lowdinOrtho(s_minushalf,F)
      # MO coeffs as linear combination of orthonormal MO basis functions
      e, Cp = eigen(Hermitian(Fp))
      # New MO Coefficients as linear combination of AO basis functions
      C = s_minushalf * Cp
      # Renew Density Matrix. 
      D, oldD = zeros(nao,nao), deepcopy(D)

      mocc = C[:,1:nocc]
      D  = mocc * mocc'
      ΔD = norm(D - oldD) # Frobenius norm

      ########################################################
      #9: COMPUTE THE NEW SCF ENERGY
      ########################################################

      oldEelec = Eelec
      Eelec = 0.0
      for i in 1:nao, j in 1:nao
         Eelec += D[i,j]*(hcore[i,j] + F[i,j])
      end

      ΔE   = Eelec - oldEelec
      Etot = Eelec + enuc

      printfmt(" {1:3d}  {2:14.8f}  {3:14.8f}  {4:14.8f} {5:14.8f}\n",
                 iteri,   Eelec,       Etot,      ΔE   ,    ΔD)
      ########################################################
      #10: TEST FOR CONVERGENCE: WHILE LOOP (go up)
      ########################################################
   end
   @info ifelse(iteri < itermax, "CONVERGED!!", "NOT CONVERGED!!")
   Etot,e,C
end

#########################################################################
# Project #4: The Second-Order Moller-Plesset
#             Perturbation Theory (MP2) Energy
#########################################################################

########################################################
#3: TRANSFORM THE TWO-ELECTRON INTEGRALS TO THE MO BASIS
########################################################

function ao2mo_noddy(nao,C,eri,mol)
   eri_mo = zeros(length(eri))
   for q in 1:nao, p in q:nao
      for r in 1:p
         lim = ifelse(p == r, q, r)
         for s in 1:lim
            pqrs = get_4idx(p,q,r,s)
            for i in 1:nao, j in 1:nao, k in 1:nao, l in 1:nao
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
   for j in 1:nao, i in j:nao
      ij = index(i,j)
      for l in 1:nao, k in l:nao
         ijkl = get_4idx(i,j,k,l)
         X[k,l] = X[l,k] = eri[ijkl]
      end

      Y = zeros(nao, nao)
      Y = C' * X
      X = zeros(nao, nao)
      X = Y * C
      for l in 1:nao, k in l:nao
         kl = index(k,l)
         tmp[kl,ij] = X[k,l]
      end
   end

   eri_mo = zeros(M*(M+1)÷2)

   for l in 1:nao, k in l:nao
      kl = index(k,l)
      X = zeros(nao,nao)
      Y = zeros(nao,nao)
      for j in 1:nao, i in j:nao
         ij = index(i,j)
         X[i,j] = X[j,i] = tmp[kl,ij]
      end
      Y = zeros(nao, nao)
      Y = C' *  X
      X = zeros(nao, nao)
      X = Y * C
      for j in 1:nao, i in j:nao
         klij = get_4idx(k,l,i,j)
         eri_mo[klij] = X[i,j]
      end
   end
   eri_mo
end

########################################################
#4: COMPUTE THE MP2 ENERGY
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
   for i in occupied, a in virtual,
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
   for i in occupied, a in virtual,
       j in occupied, b in virtual
      iajb = get_4idx(i,a,j,b)
      ibja = get_4idx(i,b,j,a)
      emp2 += eri_mo[iajb] * (2.0eri_mo[iajb]-eri_mo[ibja]) / (e[i]+e[j]-e[a]-e[b])
   end

   # In the basis of spin-molecular orbitals:
   #emp2 += 0.25eri_mo[ijab]*eri_mo[ijab]/(e[i]+e[j]-e[a]-e[b])
   emp2
end

#########################################################################
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

println( "           E(HF)     ")
println( "        -------------")
printfmt(" ME:    {1:13.10f}  \n", e_rhf_me)
printfmt(" PySCF: {1:13.10f}  \n", e_rhf_pyscf)

eri   = mol.intor("cint2e_sph", aosym="s8")
println()
e_mp2_me2   = @time mymp2(e,C,eri,mol,alg=:ao2mo_noddy) # Naive  ao2mo with my MP2
e_mp2_me3   = @time mymp2(e,C,eri,mol,alg=:ao2mo_smart) # Smart  ao2mo with my MP2
e_mp2_me1   = @time short_mp2(e,C,eri,mol)              # einsum ao2mo with my MP2
e_mp2_pyscf = @time mp.MP2(mf).kernel()[1]              # pure PySCF MP2 solution
println()

println( "           ΔE(MP2)   ")
println( "        -------------")
printfmt(" ME1:   {1:13.10f}  \n", e_mp2_me1)
printfmt(" ME2:   {1:13.10f}  \n", e_mp2_me2)
printfmt(" ME3:   {1:13.10f}  \n", e_mp2_me3)
printfmt(" PySCF: {1:13.10f}  \n", e_mp2_pyscf)

