
# Sr₂RuO₄ is a superconductor like La₍₂₋ₓ₎BaₓCuO₄ (≡ structure)
# *but* without impurities, no Cu and is ferro vs. antiferro @ T > Tc
# Intro: 10.1038/d41586-019-03734-7
# Discovery: 10.1038/372532a0
# Band calculation: 10.1103/PhysRevB.51.1385
# Fermi surface: 10.1103/PhysRevLett.85.5194
# Strain & Fermi surface: 10.1103/PhysRevLett.116.197003
# Uniaxial pressure: 10.1126/science.aaf9398
# Slater-Koster TB: 10.1103/PhysRev.94.1498
# Minimalist model with the following assumptions:
# 1. Born-Oppenheimer
# 2. Monodeterminantal wavefunction
# 3. Due to the exceptional elongation along the [001] axis
#    the system can be modeled as RuO₂ layers ⟂ [001]
#    See flat bands along Γ-X ([001] direction) @ Oguchi paper
# 4. only O 2p and Ru 4d orbitals contribute near the Fermi energy
# 5. environment simulated with a crystal field spliting
# 6. Parametrized hopping integrals (e.g. fit to ARPES)

# Default parameters
#Ep_sph = -4.41; Ed_sph = -0.41
#DQ = 0.1;       DS = -0.05;     DT = 0
#ppσ =  0.7;     ppπ = -0.3;     pdσ = -2.4; pdπ =  1

nbands = 11 # 5 Ru4d + 6 O2p orbitals
# Spherical symmetry
Ep_sph = -4.41
Ed_sph = -0.41
# Cristal field split
DQ = 0.099 # ↑ → break γ circle: Fig 1.a) 0.094 b) 0.1 c)  0.102 d) 0.103
DS = -0.05 # ↑ → break γ circle: 
DT = 0.0 # # ↑ → break γ circle
# Slater-Koster parameters
ppσ =  0.9 # ↑ → more cicular γ surface
ppπ = -0.6 # ↑ → greater γ surface
pdσ = -2.4 # no cambia nada en ∈ [2,3]
pdπ =  0.993

δEFermi = 0.01
nkmesh1d = 100 # k-mesh in each axis to integrate in the Brillouin zone

using LinearAlgebra
using Plots

include("dos_broadening.jl")

function build_H_k_RuO2layer(Ep_sph,Ed_sph,DQ,DS,DT,ppσ,ppπ,pdσ,pdπ,kx,ky)
   # Build a Hamiltonian H(kx,ky) of the set {H(kx,ky)}
   H_k = zeros(11,11)
   H_k[1,1]    = Ep_sph # O₁ 2px
   H_k[2,2]    = Ep_sph # O₁ 2py
   H_k[3,3]    = Ep_sph # O₁ 2pz
   H_k[4,4]    = Ep_sph # O₂ 2px
   H_k[5,5]    = Ep_sph # O₂ 2py
   H_k[6,6]    = Ep_sph # O₂ 2pz
   H_k[7,7]    = Ed_sph + 6DQ - 2DS - 6DT # Ru 4dz²    |Ru atoms lie in a distorted octahedron
   H_k[8,8]    = Ed_sph + 6DQ + 2DS -  DT # Ru 4dx²-y² |thus, d → t_2g ⊕ e_g with DQ gap
   H_k[9,9]    = Ed_sph - 4DQ + 2DS -  DT # Ru 4dxy    |There is also a tetragonal distortion
   H_k[10,10]  = Ed_sph - 4DQ -  DS + 4DT # Ru 4dyz    |induced by uniaxial strain, gaps DT and DS
   H_k[11,11]  = Ed_sph - 4DQ -  DS + 4DT # Ru 4dzx    |Splitings maintain the baricenter.
   # See Table I. in Slater-Koster article.
   # NOTE: trick of p -> ip to have a real matrix instead of complex
   # Top view ⟂ [001]:
   #   Ru   O1   Ru   O1   Ru
   #           + - ------+
   #   O2      | O2      | O2
   #           |         | 
   #   Ru   O1 | Ru   O1 | Ru
   #           +---------+
   #   O2        O2        O2
   #
   #   Ru   O1   Ru   O1   Ru
   H_k[4,1]    =  2(ppσ+ppπ)*cos(0.5kx)*cos(0.5ky)
   H_k[5,1]    = -2(ppσ-ppπ)*sin(0.5kx)*sin(0.5ky)
   H_k[4,2]    = -2(ppσ-ppπ)*sin(0.5kx)*sin(0.5ky)
   H_k[5,2]    =  2(ppσ+ppπ)*cos(0.5kx)*cos(0.5ky)
   H_k[6,3]    =  4ppπ*cos(0.5kx)*cos(0.5ky)
   H_k[7,1]    = -pdσ*sin(0.5kx)
   H_k[8,1]    =  √3*pdσ*sin(0.5kx)
   H_k[9,2]    =  2pdπ*sin(0.5kx)
   H_k[11,3]   = -2pdπ*sin(0.5kx)
   H_k[9,4]    =  2pdπ*sin(0.5ky)
   H_k[7,5]    = -pdσ*sin(0.5ky)
   H_k[8,5]    = -√3*pdσ*sin(0.5ky)
   H_k[10,6]   = -2pdπ*sin(0.5ky)
   H_k = Symmetric(H_k,:L)
end

# ----------------------------------------------------------------------------------------
# Band: k-path
kpath = zeros(2,342) # {(kx ky)ᵀ}
kpath[1,1:101]    = range(0,stop=π,length=101)
kpath[2,1:101]   .= 0

kpath[1,102:201] .= π
kpath[2,102:201]  = range(0,stop=π,length=100)

kpath[1,202:342]  = range(π,stop=0,length=141)
kpath[2,202:342]  = range(π,stop=0,length=141)

nkpoints = size(kpath)[2]
Enk = zeros(nbands,nkpoints)

for k in 1:nkpoints
   kx = kpath[1,k]
   ky = kpath[2,k]
   H_k = buildH_RuO2layer(Ep_sph,Ed_sph,DQ,DS,DT,ppσ,ppπ,pdσ,pdπ,kx,ky)
   εk, ψk = eigen(H_k)
   Enk[:,k] = εk
end

emin, emax = extrema(vec(Enk))
band=plot(1:nkpoints,[Enk[n,:] for n in 1:nbands],
          xticks=([1,101,201],["(0,0)", "(1,0)", "(1,1)"]),
          ylims=(emin,emax),
          #label=permutedims(["Band $n" for n in 1:nbands]), # must be row vector
          leg=false)


# ----------------------------------------------------------------------------------------
# DOS, PDOS & Fermi surface: Brillouin zone integration
nkpoints = nkmesh1d*nkmesh1d
Enk = zeros(nbands,nkmesh1d*nkmesh1d)
PDOS = zeros(nkmesh1d*nkmesh1d,nbands,nbands)

kx_Fermi_surf = []
ky_Fermi_surf = []

for nx in 1:nkmesh1d, ny in 1:nkmesh1d
   kx = π*(nx-0.5)/nkmesh1d
   ky = π*(ny-0.5)/nkmesh1d
   H_k = buildH_RuO2layer(Ep_sph,Ed_sph,DQ,DS,DT,ppσ,ppπ,pdσ,pdπ,kx,ky)
   εk, ψk = eigen(H_k)
   if minimum(abs.(εk)) < δEFermi
      push!(kx_Fermi_surf,kx/π)
      push!(ky_Fermi_surf,ky/π)
   end
   kindex = (nx-1)*nkmesh1d + ny
   Enk[:,kindex] = εk 
   for n in 1:nbands
      PDOS[kindex,n,:] = ψk[:,n].*ψk[:,n] # weight of AO to state
   end
end

dos, e_dos = dos_broadening(vec(Enk))
dos = plot(e_dos, dos,xlabel="DOS",leg=false,
           ylims=(emin,emax))
l = @layout [a b]
plot(band,dos,layout=l)
savefig("Sr2RuO4_banddos.pdf")

pdos = plot(title="Sr2RuO4 PDOS",xlabel="Energy")
for n in 1:nbands
   plot!(pdos,dos_broadening(vec(Enk),vec(PDOS[:,:,n]))...,
        label="AO $n")
end
savefig(pdos,"Sr2RuO4_pdos.pdf")

fermi_surf = scatter(title="Sr2RuO4 Fermi surface",leg=false,
             xticks=([-1,-0.5,0,0.5,1],["-\\pi", "-\\pi/2", "0", "\\pi/2", "\\pi"]),
             yticks=([-1,-0.5,0,0.5,1],["-\\pi", "-\\pi/2", "0", "\\pi/2", "\\pi"]),
             xlabel = "k_x", ylabel = "k_y", size=(400,400)
            )
scatter!(fermi_surf, kx_Fermi_surf, ky_Fermi_surf,markercolor=:red,markersize=1.7,markerstrokewidth=0)
scatter!(fermi_surf,-kx_Fermi_surf, ky_Fermi_surf,markercolor=:red,markersize=1.7,markerstrokewidth=0)
scatter!(fermi_surf, kx_Fermi_surf,-ky_Fermi_surf,markercolor=:red,markersize=1.7,markerstrokewidth=0)
scatter!(fermi_surf,-kx_Fermi_surf,-ky_Fermi_surf,markercolor=:red,markersize=1.7,markerstrokewidth=0)
savefig(fermi_surf,"Sr2RuO4_Fermi_surface.pdf")

