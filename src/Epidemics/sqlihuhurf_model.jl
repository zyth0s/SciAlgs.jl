
# Modelización Epidemiológica del Covid-19 para España 1 de Abril de 2020
# Simón, Clara Burgos
# Cortés, Juan-Carlos
# Navarro, Elena López
# Martínez-Rodríguez, David
# Martínez-Rodríguez, Pablo
# Julián, Raul S.
# Villanueva, Rafael-J.
# Instituto Universitario de Matemática Multidisciplinar,
# Universitat Politècnica de València, 46022 Valencia, España.
#
# Empieza el 27 de enero (incluido)

using DifferentialEquations
using Plots

# Model
sqlihuhurf_ode = @ode_def SQLIHUHURFModel begin
  dS  = -β(t)*S*I/Pop - δ(t)*S + τ*Q
  dQ  = δ(t)*S - τ*Q
  dL  = β(t)*S*I/Pop - γ₁*L
  dI  = γ₁*L - (γ₂ + α₁)*I
  dH  = γ₂*I - (d₁ + α₂ + γ₃)*H
  dU  = γ₃*H - (d₂ + α₃)*U
  dHU = α₃*U - η*HU
  dR  = α₁*I + α₂*H + η*HU
  dF  = d₁*H + d₂*U
end Pop τ γ₁ γ₂ γ₃ α₁ α₂ α₃ η d₁ d₂

β(t) = t <= 17 ? 0.59   : 0.1 # 14 de marzo se declara el estado de alarma
δ(t) = t == 19 ? 0.7036 : 0.0 # 16 de marzo se hace efectivo el estado de alarma

# Parameters
Pop = float(47100396) # Total population in Spain
#β₀  = 0.59 # 0.1 despues del estado de alarma
#δ   = 0.7036
p₁  = 0.05
p₂  = 0.065
p₃  = 0.045
p₄  = 0.55
L₀  = 10926.0
I₀  = 15194.0
#
τ = 0
γ₁ = 1/5.2
γ₂ = p₁/5.8
γ₃ = p₃
α₁ = (1-p₁)/14
α₂ = (1-p₂-p₃)/7
α₃ = (1-p₄)/14
η = 1.0/6
d₁ = p₂/7.5
d₂ = p₄/8


par=[Pop, τ, γ₁, γ₂, γ₃, α₁, α₂, α₃, η, d₁, d₂]
init=[Pop-L₀-I₀, 0, L₀, I₀, 0, 0, 0, 0, 0]
tspan=(0.0,61)

sqlihuhurf_prob = ODEProblem(sqlihuhurf_ode,init,tspan,par)
#sol = solve(sqlihuhurf_prob,saveat=1.0,reltol=1e-8,abstol=1e-8,alg_hints=[:stiff])
sol = solve(sqlihuhurf_prob,saveat=1.0,reltol=1e-8,abstol=1e-8)


#plot(sol)
#plot(sol.t,sol[1,:],label="susceptibles")
plot(sol.t,sol[4,:],label="infected")
#plot(sol.t,sol[5,:],label="hospitalized")
#plot(sol.t,sol[9,:],label="deaths")
#plot(sol.t,[sol[5,:],sol[4,:]])

