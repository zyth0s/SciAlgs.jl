
# https://simonverret.github.io/2018/11/15/visual-greens-functions.html
# Two level system with Green's functions

using LinearAlgebra: SymTridiagonal,eigvals, eigvecs, Diagonal, diag


Δ = 0.2  # The gap
η = 1e-2 # infinitesimal
Nk = 100
Nϵ = 100

ϵmin, ϵmax = -2,2
kmesh = range(-π,stop=π,length=Nk)
ϵmesh = range(ϵmin,stop=ϵmax,length=Nϵ)

ε(k) = [k,-k] # crossing levels

H(k) = SymTridiagonal(ε(k),[Δ])

G(k,ϵ) = inv(Diagonal(fill(ϵ,2)) - H(k) + Diagonal(fill(η*im,2)))

# In the diagonal basis we can select the Green function for each band E^±
# G^±(k,ϵ) = inv(ϵ - ε(k) + iη)
function Gdiag(k,ϵ) 
   Epm = eigvals(H(k))
   inv(Diagonal(fill(ϵ,2)) - Diagonal(Epm) + Diagonal(fill(η*im,2)))
end

# The spectral weight A(k,ϵ) = -Im{G(k,ϵ)}
# equivalently, A^±(k,ϵ) = -Im{ inv(ϵ - ε(k) + iη) }
A(k,ϵ)     = -imag( G(k,ϵ)     )
Adiag(k,ϵ) = -imag( Gdiag(k,ϵ) )

# The density of states
# N(ϵ)   = 1/π ∫ -Im{G(k,ϵ)}   dk = 1/π ∫ A(k,ϵ)   dk
# N^±(ϵ) = 1/π ∫ -Im{G^±(k,ϵ)} dk = 1/π ∫ A^±(k,ϵ) dk

using PyPlot

# Energy dispersion curves of two levels
E = eigvals.(H.(kmesh))

subplot(121) # 1rows 2cols 1st plot
xlabel(L"k")
ylabel(L"E_k^\pm")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
#axvline(color="k",ls="--",lw=0.5)
#axhline(color="k",ls="--",lw=0.5)
plot(kmesh,[elem[1] for elem in E],label=L"E^- Ground")
plot(kmesh,[elem[2] for elem in E],label=L"E^+ Excitd")
plot(kmesh,ε.(kmesh),lw=0.1,color="k", label=[L"ϵ_1", L"ϵ_2"])
#legend(loc="lower center")

subplot(122) # 1rows 2cols 2nd plot
# N(ϵ) = 1/π ∫ A(k,ϵ) dk
_kmesh = range(-π,stop=π,length=6Nk) # NOTE: fine grid needed
pdos11 = (1/π)*sum([A(k,ϵ)[1,1] for ϵ in ϵmesh, k in _kmesh],dims=2)/Nk * 2
#pdos22 = (1/π)*sum([A(k,ϵ)[2,2] for ϵ in ϵmesh, k in _kmesh],dims=2)/Nk * 2
xlabel(L"N(ϵ)")
ylabel(L"ϵ")
plot(pdos11,ϵmesh,label="pdos11",color="k")
#plot(pdos22,ϵmesh,label="pdos22",color="k")
# N^±(ϵ) = 1/π ∫ A^±(k,ϵ) dk
pdos1  = (1/π)*sum([Adiag(k,ϵ)[1,1] for ϵ in ϵmesh, k in _kmesh],dims=2)/Nk
pdos2  = (1/π)*sum([Adiag(k,ϵ)[2,2] for ϵ in ϵmesh, k in _kmesh],dims=2)/Nk
fill_betweenx(ϵmesh,vec(pdos1),label="pdos-")
fill_betweenx(ϵmesh,vec(pdos2),label="pdos+")
tight_layout()
savefig("../figures/two_level_bands_dos.pdf")
clf()

# Plot eigenvectors
V = eigvecs.(H.(kmesh))

subplot(211) # 2rows 2cols 1st plot
ylabel(L"v_1^\pm")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
axvline(color="k",ls="--",lw=0.5)
axhline(color="k",ls="--",lw=0.5)
plot(kmesh,[v[1,1] for v in V],label="v1-")
plot(kmesh,[v[2,1] for v in V],label="v1+")

subplot(212) # 2rows 2cols 2nd plot
xlabel(L"k")
ylabel(L"v_2^\pm")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
axvline(color="k",ls="--",lw=0.5)
axhline(color="k",ls="--",lw=0.5)
plot(kmesh,[v[1,2] for v in V],label="v2-")
plot(kmesh,[v[2,2] for v in V],label="v2+")
# FIXME: incorrect curves
savefig("../figures/two_level_eigvecs.pdf")
clf()

# Colormap of Im{ G(k,ϵ)[1,1] }
subplot(221)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [imag(G(k,ϵ)[1,1]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Im{ G(k,ϵ)[2,1] }
subplot(222)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [imag(G(k,ϵ)[2,1]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Im{ G(k,ϵ)[1,2] }
subplot(223)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [imag(G(k,ϵ)[1,2]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Im{ G(k,ϵ)[2,2] }
subplot(224)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [imag(G(k,ϵ)[2,2]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)
# Colorbar
plt.subplots_adjust(bottom=0.1,right=0.8,top=0.9)
cax = plt.axes([0.85,0.1,0.075,0.8])
colorbar(cax=cax,ticks=collect(-1:1:1))
#tight_layout()
savefig("../figures/two_level_Gimag.pdf")
clf()


# Colormap of Re{ G(k,ϵ)[1,1] }
subplot(221)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [real(G(k,ϵ)[1,1]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Re{ G(k,ϵ)[2,1] }
subplot(222)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [real(G(k,ϵ)[2,1]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Re{ G(k,ϵ)[1,2] }
subplot(223)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [real(G(k,ϵ)[1,2]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)

# Colormap of Re{ G(k,ϵ)[2,2] }
subplot(224)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Grealtmp = [real(G(k,ϵ)[2,2]) for ϵ in ϵmesh, k in kmesh]
imshow(Grealtmp,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower",vmin=-1,vmax=1)
# Colorbar
plt.subplots_adjust(bottom=0.1,right=0.8,top=0.9)
cax = plt.axes([0.85,0.1,0.075,0.8])
colorbar(cax=cax,ticks=collect(-1:1:1))
#tight_layout()
savefig("../figures/two_level_Greal.pdf")

