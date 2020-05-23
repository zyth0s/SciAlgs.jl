
# Adapted from https://hub.gke.mybinder.org/user/marcelloperathoner-covid-19-wqgzqk30/notebooks/jupyter/seirc.ipynb
using DifferentialEquations
using Plots

seirc_ode = @ode_def SERICModel begin
    F = β * I / N
    dS = -F * S
    dE = F * S - α * E
    dI = α * E - γ * I
    dR = γ * I
    dC = c1 * α * E - c2 * C
    end α β γ c1 c2 N

# COVID-19 data
N = 83019213.0      # Germany: Population 31.12.2018
S = 4838.0          # Germany: Confirmed Covid-19 cases 16.03.2020

α = 2    # How many days before E(xposed) becomes I(nfectious)
β = 0.2  # How many people does an I infect each day on average
γ = 7    # How many days an I stays infectious on average

c1 = 0.0165 # How many I become C(ritical)
c2 = 28     # How many days a C occupies a bed
c3 = 5600   # How many beds thare are (28000 * 20% free)

oversample = 10 # how many iterations per simulated day
parms = [1 / α, β, 1 / γ, c1, 1 / c2, N] 
initial  = [N - S, S, 0, 0, 0]     # S, E, I, R, C
times = range(0, stop=365, length=365 * oversample + 1)
tspan = (0.0, 365) 


seirc_prob = ODEProblem(seirc_ode,initial,tspan,parms)
seirc_sol = solve(seirc_prob,saveat = 1.0/oversample );

plot(seirc_sol,xlabel="Time",ylabel="Number")
plot!(seirc_sol.t,ones(365*oversample+1)*c3, label="Beds")

# Magnify the peak demand of hospital beds
ylims!(0,c3*1e2)
