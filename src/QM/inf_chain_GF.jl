# Julia programming language version 1.0.2
# Odashima, Prado, Vernek. "Pedagogical introduction to equilibrium Green's functions: 
#                           condensed-matter examples with numerical implementations"
#                           10.1590/1806-9126-RBEF-2016-0087
# §2.4 Infinite linear chain (solved with Green's functions method)

using PyPlot                        # Matplotlib library

ϵ₀ = 0                              # local site energy
η = eps()                           # positive infinitesimal
ωmin = -2; ωmax = 2                 # energy range
Nω = 1000                           # number of energy points
ω = range(ωmin,stop=ωmax,length=Nω) # vector of energies
g = @. 1 / (ω - ϵ₀ + η*im)          # undressed propagator, eq. (66)
t = 1                               # symmetric real hopping

# Semi-infinite chain analytic expression G₁₁, eq. (78)

Gsemi = @. (1 / (2g*t^2)) * (1 - sqrt(1-4t^2*g^2))

# Infinite chain analytic expression obtained
# by joining two semi-infinite chains, eq. (83)

Ginf = @. Gsemi/(1-Gsemi^2*t^2)

xlabel(L"Energy, ($\omega - \epsilon_0 $)/t")
ylabel(L"Density of states, $\rho_{11} (\omega) $")
axis([ωmin,ωmax,0,1.4])
plot(ω, (-1.0/π)*imag(Ginf), linewidth=3.0) # eq. (44)
