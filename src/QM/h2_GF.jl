# Julia programming language version 1.0.2
# Odashima, Prado, Vernek. "Pedagogical introduction to equilibrium Green's functions: 
#                           condensed-matter examples with numerical implementations"
#                           10.1590/1806-9126-RBEF-2016-0087
# §2.2 Two-site chain: H₂ with tight-binding hamiltonian
# solved with retarded Green's functions method

using PyPlot                        # Matplotlib library

ϵ₀ = 0                              # local site energy
η = 0.01                            # positive infinitesimal
ωmin, ωmax = -2, 2                  # energy range
Nω = 1000                           # number of energy points
ω = range(ωmin,stop=ωmax,length=Nω) # vector of energies
g = @. 1 / (ω - ϵ₀ + η*im)          # undressed propagator, eq. (66)
t = 1                               # symmetric real hopping

# H₂ analytic expression for G₁₁, eq. (72)

G11 = @. g / (1 - t^2*g^2)

# Figure 2 in paper
xlabel(L"Energy, ($\omega - \epsilon_0 $)/t")
ylabel(L"Local density of states, $\rho_{11} (\omega) $")
axis([ωmin,ωmax,0,16])
fill(ω, (-1.0/π)*imag(G11))         # eq. (44)
