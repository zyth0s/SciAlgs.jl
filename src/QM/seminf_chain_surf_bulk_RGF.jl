# Julia programming language version 1.0.2
# Odashima, Prado, Vernek. "Pedagogical introduction to equilibrium Green's functions: 
#                           condensed-matter examples with numerical implementations"
#                           10.1590/1806-9126-RBEF-2016-0087
# §3.2 Semi-infinite linear chain (with surface-bulk recursive Green's functions)

using PyPlot                        # Matplotlib library

function semi_inf_chain_surf_bulk_recursive(Ndec=16, # number of decimation iterations
                                            η=eps()) # positive infinitesimal
   ϵ₀ = 0                              # local site energy
   ωmin, ωmax = -2, 2                  # energy range
   Nω = 1000                           # number of energy points
   ω = range(ωmin,stop=ωmax,length=Nω) # vector of energies
   g = @. 1 / (ω - ϵ₀ + η*im)          # undressed propagator, eq. (66)
   g10 = g20 = g30 = g                 # initialization of undressed retarded GF
   t = td = ones(Nω)                   # symmetric real hopping

   for i in 1:Ndec # Decimation Loop

      g1 = @. g10/(1 - g10*t*g20*td) # effective Green’s function of site 1, eq. (123)
      g2 = @. g20/(1 - g20*td*g20*t - g20*t*g20*td) # effective Green’s function of site 2, eq. (129)
      g3 = @. g30/(1 - g30*td*g20*t) # eq. (132)

      t  = t  .* g20 .* t            # Renormalization of the hoppings, t̃
      td = td .* g20 .* td           # Note that we do not conjugate g20, t̃*

      g10 = g1                       # Update of the loop variables
      g20 = g2
      g30 = g3

   end
   G11 = @. g10/(1 - g10*t*g20*td/(1 - g20*t*g30*td)) # final surface Green’s function of site 1, eq. (124)

   ω, G11
end

# Figure 13 in paper
xlabel(L"Energy, ($\omega - \epsilon_0 $)/t")
ylabel(L"Local density of states, $\rho_{11} (\omega) $")
axis([-2,2,0,8])
for Ndec in 0:3
   ω, G11 = semi_inf_chain_surf_bulk_recursive(Ndec,0.02)

   plot(ω, (-1/π)*imag(G11), label="$(Ndec)th iteration: $(2^(Ndec+1)+1) sites")  # Plotting the density of states of the surface site, eq. (44)
end
legend()

clf()
# Figure 14 in paper
title("Semi-infinite chain with surface-bulk RFG")
xlabel(L"Energy, ($\omega - \epsilon_0 $)/t")
ylabel(L"Local density of states, $\rho_{11} (\omega) $")
axis([-2,2,0,1.6])
for Ndec in 11:14
   ω, G11 = semi_inf_chain_surf_bulk_recursive(Ndec,1e-4)

   plot(ω, (-1/π)*imag(G11), label="$(Ndec)th iteration: $(2^(Ndec+1)+1) sites")  # Plotting the density of states of the surface site, eq. (44)
end
legend()
