
using PyCall: pyimport
using Formatting: printfmt
using LinearAlgebra: eigen, Diagonal, dot, norm, Hermitian

pyscf = pyimport("pyscf")

mol = pyscf.gto.Mole()
mol.build(verbose = 0,
          #atom = join(split(read(open("../data/h2.xyz"),String),"\n")[3:end],"\n"), # from XYZ file format
          #atom = join(split(read(open("../data/acetaldehyde.xyz"),String),"\n")[3:end],"\n"), # from XYZ file format
          atom = """8 0 0. 0
          1 0 -0.757 0.587
          1 0 0.757 0.587""",
          basis = "sto-3g",
         )

index(i,j) = max(i-1,j-1)*(max(i-1,j-1)+1)÷2 + min(i-1,j-1) + 1 # 1 based indexing

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
    F = zeros(nao,nao)  
    for i in 1:nao, j in 1:nao
        F[i,j] = hcore[i,j]  
        for k in 1:nao, l in 1:nao
            ijkl = get_4idx(i,j,k,l)
            ikjl = get_4idx(i,k,j,l)
            F[i,j] += D[k,l]*(2eri[ijkl]-eri[ikjl])  
         end
    end
    F
end


function scf_rhf(mol)
   nelec = mol.nelectron
   nocc  = mol.nelectron ÷ 2 
   nao   = mol.nao_nr() |> n -> convert(Int,n) # why not converted ??!
   M     = nao*(nao+1)÷2
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
                 iteri,   Eelec,       Etot,        ΔE   ,     ΔD)
      ########################################################
      #10: TEST FOR CONVERGENCE: WHILE LOOP (go up)
      ########################################################
   end
   @info ifelse(iteri < itermax, "CONVERGED!!", "NOT CONVERGED!!")
   Etot
end

e_rhf_me    =       scf_rhf(mol)
e_rhf_pyscf = pyscf.scf.RHF(mol).kernel()

println( "           E(HF)     ")
println( "        -------------")
printfmt(" ME:    {1:13.10f}  \n", e_rhf_me)
printfmt(" PySCF: {1:13.10f}  \n", e_rhf_pyscf)
