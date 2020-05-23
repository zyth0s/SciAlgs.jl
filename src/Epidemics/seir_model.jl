# http://epirecip.es/epicookbook/chapters/seir/julia
# SEIR model
# S is the fraction of susceptible individuals (able to contract the disease)
# E is the fraction of exposed individuals (infected but not yet infectious)
# I is the fraction of infective individuals (capable of transmitting the disease)
# R is the fraction of recovered individuals (became permanently immune)
# S + E + I + R = 1
# μ is the birth/dead rate
# 1/α is the mean latent period for the disease
# 1/γ is the mean infectious period
# β is the contact rate, assumed constant
# R₀ = βα / (μ+α)(μ+γ) is the reproductive number
#
using DifferentialEquations
using Plots

seir_ode = @ode_def SEIRModel begin
  dS = μ - β*S*I - μ*S
  dE = β*S*I - (α+μ)*E
  dI = α*E - (γ+μ)*I
  end β α γ μ

par=[520/365,1/60,1/30,774835/(65640000*365)]
init=[0.99,0.005,0.005]
tspan=(0.0,365.0)

seir_prob = ODEProblem(seir_ode,init,tspan,par)

sol = solve(seir_prob)

# R = 1 - S - E - I
R=ones(1,size(sol,2))-sum(sol,dims=1)

plot(sol.t,[sol',R'],xlabel="Time",ylabel="Proportion", label = ["S","E","I","R"])

