
# https://simonverret.github.io/2018/11/15/visual-greens-functions.html
# Free electron system with Green's functions

# atomic units
hbar = 1
m = 1

Nk = 100
Nϵ = 100
kmesh = range(-π,stop=π,length=Nk)
η = 1e-1
ϵ = 0.0

ε(k) = hbar^2 * k^2 / (2m)

G(k,ϵ) = 1/(ϵ - ε(k) + η*im)

using PyPlot

# Plot of energy dispersion ε(k)
#xlabel(L"k")
#ylabel(L"\varepsilon (k)")
#plot(kmesh, ε.(kmesh))
#grid("on")

ϵmesh = range(-1,stop=maximum(ε.(kmesh)),length=Nϵ)


# Colormap of Re{ G(k,ϵ) }
subplot(221)
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Greal = real.([G(k,ϵ) for ϵ in ϵmesh, k in kmesh])
imshow(Greal,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower")
colorbar()
annotate(xycoords="axes fraction",xy=(0.4,0.5),L"Re $G$")
clim([-1,1])

# Colormap of -Im{ G(k,ϵ) }
subplot(223)
xlabel(L"k")
ylabel(L"ϵ")
xticks([-π,-π/2,0,π/2,π],["-π","-π/2","0","π/2","π"])
Gimag = imag.([G(k,ϵ) for ϵ in ϵmesh, k in kmesh])
imshow(-Gimag,cmap="RdBu",extent=[kmesh[1],kmesh[end],ϵmesh[1],ϵmesh[end]],
       origin="lower")
colorbar()
annotate(xycoords="axes fraction",xy=(0.4,0.5),L"-Im $G$")
clim([-1,1])


# Lineplot of Re{G(k=π/2,ϵ)}
subplot(222)
ylabel(L"ϵ")
xlabel(L"Re $G(ϵ)$")
Greal = real.([G(π/2,ϵ) for ϵ in ϵmesh])
plot(Greal, ϵmesh)
axvline(color="k",ls="--",lw=0.5)

# Lineplot of -Im{G(k=π/2,ϵ)} = spectral weight
subplot(224)
ylabel(L"ϵ")
xlabel(L"-Im $G(ϵ)$")
Gimag = imag.([G(π/2,ϵ) for ϵ in ϵmesh])
plot(-Gimag, ϵmesh)
tight_layout()
savefig("../figures/free_elec_band_GF.pdf")


