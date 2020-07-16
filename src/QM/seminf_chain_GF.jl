# Julia programming language version 1.0.2
# Odashima, Prado, Vernek. "Pedagogical introduction to equilibrium Green's functions: 
#                           condensed-matter examples with numerical implementations"
#                           10.1590/1806-9126-RBEF-2016-0087
# §3.2 Semi-infinite linear chain (with surface-bulk retarded Green's functions)

using PyPlot                        # Matplotlib library


η=eps()                             # positive infinitesimal
ϵ₀ = 0                              # local site energy
ωmin, ωmax = -4, 4                  # energy range
Nω = 1000                           # number of energy points
ω = range(ωmin,stop=ωmax,length=Nω) # vector of energies
g = @. 1 / (ω - ϵ₀ + η*im)          # undressed propagator, eq. (66)
t = 1                               # symmetric real hopping

# Semi-infinite chain analytic, eq. (78), choose sign -

G11 = @. 1/(2g*t^2) * (1 - sqrt(1-4t^2*g^2))

# Figure 5 in paper; Fig 5.10 Economou 3ed
xlabel(L"Energy, ($\omega - \epsilon_0 $)/t")
axis([-4,4,-1,1])
plot(ω, -abs(t)*imag(G11), label=L"-|t| \; Im(G_{11})")
fill(ω, -abs(t)*imag(G11))
plot(ω,  abs(t)*real(G11), label=L"\quad |t| \; Re(G_{11})",ls="--")
legend(loc="lower right")

