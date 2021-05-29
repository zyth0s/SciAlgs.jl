
# Tight-binding calculations of the valence bands of diamond and zincblende crystals
# D. J. Chadi, M. L. Cohen. doi:10.1002/pssb.2220680140

using LinearAlgebra

# eqs 7-10
function compute_g(k)
   @assert size(k) == (3,)
   # Si (diamond structure)
   [cospi(0.5k[1])*cospi(0.5k[2])*cospi(0.5k[3]) - im*sinpi(0.5k[1])*sinpi(0.5k[2])*sinpi(0.5k[3]),
   -cospi(0.5k[1])*sinpi(0.5k[2])*sinpi(0.5k[3]) + im*sinpi(0.5k[1])*cospi(0.5k[2])*cospi(0.5k[3]),
   -sinpi(0.5k[1])*cospi(0.5k[2])*sinpi(0.5k[3]) + im*cospi(0.5k[1])*sinpi(0.5k[2])*cospi(0.5k[3]),
   -sinpi(0.5k[1])*sinpi(0.5k[2])*cospi(0.5k[3]) + im*cospi(0.5k[1])*cospi(0.5k[2])*sinpi(0.5k[3])]
end

# eq 6.
function build_H(g, skparams)

   Es, Ep, Vss, Vsp, Vxx, Vxy = skparams
   # Es = Es0 = Es1
   # Ep = Ep0 = Ep1
   # Vsp = Vs0p = Vs1p
   𝟎 = 1 # 1- to 0-based indexing
   𝟏 = 2
   𝟐 = 3
   𝟑 = 4
   # :AlignCtrl =Clp1P1IW \S\+
   # :Align
   #          s₀       s₁        x₀        y₀       z₀      x₁      y₁       z₁ 
   [Es        Vss*g[𝟎]   0           0           0          Vsp*g[𝟏]  Vsp*g[𝟐]  Vsp*g[𝟑];  # s₀ 
   Vss*g[𝟎]'  Es         -Vsp*g[𝟏]'  -Vsp*g[𝟐]'  -Vsp*g[𝟑]' 0         0         0;         # s₁ 
   0          -Vsp*g[𝟏]  Ep          0           0          Vxx*g[𝟎]  Vxy*g[𝟑]  Vxy*g[𝟏];  # x₀ 
   0          -Vsp*g[𝟐]  0           Ep          0          Vxy*g[𝟑]  Vxx*g[𝟎]  Vxy*g[𝟏];  # y₀ 
   0          -Vsp*g[𝟑]  0           0           Ep         Vxy*g[𝟏]  Vxy*g[𝟐]  Vxx*g[𝟎];  # z₀ 
   Vsp*g[𝟏]'  0          Vxx*g[𝟎]'   Vxy*g[𝟑]'   Vxy*g[𝟏]'  Ep        0         0;         # x₁ 
   Vsp*g[𝟐]'  0          Vxy*g[𝟑]'   Vxx*g[𝟎]'   Vxy*g[𝟐]'  0         Ep        0;         # y₁ 
   Vsp*g[𝟑]'  0          Vxy*g[𝟏]'   Vxy*g[𝟏]'   Vxx*g[𝟎]'  0         0         Ep]        # z₁ 
   #   |
   #   g[𝟑]?
end

function eigenvalues(g,skparams)
   H = build_H(g,skparams)
   @assert H ≈ Hermitian(H)
   eigvals(Hermitian(H))
end

# eq. 12
function energy_s_decoupled(k, skparams)
   Es, _, Vss, = skparams
   g0 = compute_g(k)[1]
   Es +  Vss*abs(g0), Es - Vss*abs(g0)
end

function bands(kpath, skparams)
   bands = []
   for k in kpath
      g = compute_g(k)
      e = eigenvalues(g, skparams)
      push!(bands, e)
   end
   bands
end

function build_kpath(k0, k1, length=100)
   k01 = k1-k0 # vector pointing from k0 to k1
   kpath = []
   for fraction in range(0,1,length=length)
      push!(kpath, k0 + fraction * k01)
   end
   kpath
end

# Example: Si (diamond structure)

# Lattice parameter
a = 5.43071 # Å; ICSD:29287 doi:10.1016/0022-3697(60)90069-X

L = 2π*[0.5, 0.5, 0.5]/a
Γ = 2π*[0.0, 0.0, 0.0]/a
X = 2π*[0.0, 0.0, 1.0]/a
# U
# K
# T

LΓX = vcat(build_kpath(L,Γ),
           build_kpath(Γ,X))

# from Table 3 (no 2ⁿᵈ nearest neighbor interaction)
skparams = (Es = 0, Ep = 7.20, Vss = -8.13, Vsp = 5.88, Vxx = 3.17, Vxy = 7.51)

_bands = bands(LΓX, skparams)

using PyPlot
for iband in 1:length(_bands[1])
   plot(getindex.(_bands,iband), label="Band $iband")
end
# Approximation of lowest state
e_lowest_approx = [energy_s_decoupled(k,skparams)[1] for k in LΓX]
#e_second_approx = [energy_s_decoupled(k,skparams)[2] for k in LΓX]
plot(e_lowest_approx, "--", label="Approx. lowest")
#plot(e_second_approx, "--", label="Approx. second")

legend() #bbox_to_anchor=(1,0), loc="lower left")
title("Si")
ylabel("Energy (eV)")
xlabel("k")
xticks([1,100,200], ["L", "Γ", "X"])
savefig("../figures/silicon_bands_chadi_cohen.pdf")

E_valence_max    = maximum(getindex.(_bands,4))
E_conduction_min = minimum(getindex.(_bands,5))

@info "" E_valence_max
@info "" E_conduction_min
