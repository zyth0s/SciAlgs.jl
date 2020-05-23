
# http://epirecip.es/epicookbook/chapters/sis/intro
using DifferentialEquations
using Plots

sis_ode = @ode_def SISModel begin
    dS = -β*S*I + γ*I
    dI =  β*S*I - γ*I
    end β γ

#parms = [0.1,0.05] # β, γ
#
# COVID-19
# incubation + delay to hospitalization + avg. hospitalization.
D = 5.1 + 5 + 10.4 # average days to recover = 1/γ;
R₀ = 2.2 # = β/γ
#--
γ = 1/D
β = R₀ * γ
HIT = 1 - 1/R₀ # herd immune threshold
parms = [β, γ]
init = [0.99,0.01]
tspan = (0.0,200.0)

sis_prob = ODEProblem(sis_ode,init,tspan,parms)
sis_sol = solve(sis_prob,saveat = 0.1);


plot(sis_sol,xlabel="Time",ylabel="Fraction")
