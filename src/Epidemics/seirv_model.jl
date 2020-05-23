
# A mathematical model for the novel coronavirus epidemic in Wuhan, China; Chayu Yang and Jin Wang∗
# doi: 10.3934/mbe.2020148
# Discussion:
# https://francis.naukas.com/2020/03/14/el-modelo-seirv-aplicado-a-la-epidemia-de-coronavirus-en-wuhan-china/
using DifferentialEquations
using Plots

# Model
seirv_ode = @ode_def SEIRVModel begin
   dS = Λ - βE(E)*S*E - βI(I)*S*I - βV(V)*S*V - μ*S
   dE = βE(E)*S*E + βI(I)*S*I + βV(V)*S*V - (α + μ)*E
   dI = α*E - (w + γ + μ)*I
   dR = γ*I - μ*R
   dV = ξ₁*E + ξ₂*I - σ*V
end α γ μ Λ w σ ξ₁ ξ₂

βE(E) = βE₀/(1+c*E)
βI(I) = βI₀/(1+c*I)
βV(V) = βV₀/(1+c*V)

# Model parameters
S = 8998505
E = 1000
I = 475
R = 10
V = 10000.0
Λ   = 271.23   # Influx rate
βE₀ = 3.11e-8  # Transmission constant between S and E
βI₀ = 0.62e-8  # Transmission constant between S and I
βV₀ = 1.03e-8  # Transmission constant between S and V; fitted
c   = 1.01e-4  # Transmission adjustment coefficient; fitted
μ   = 3.01e-5  # Natural death rate
α   = 1.0/7.0  # 1/Incubation period
w   = 0.01     # Disease-induced death rate
γ   = 1.0/15.0 # Recovery rate
σ   = 1.0      # Removal rate of virus
ξ₁  = 2.30     # Virus shedding rate by exposed people; fitted
ξ₂  = 0.0      # Virus shedding rate by infected people

params=[α, γ, μ, Λ, w, σ, ξ₁, ξ₂]
initial=[S, E, I, R, V]
tspan=(0.0,365.0*8.21)

# Epidemic curve prediction
seirv_prob = ODEProblem(seirv_ode,initial,tspan,params)
sol = solve(seirv_prob,alg=BS3(),saveat=1.0)

l = @layout [a b]
p1 = plot(sol.t[1:300],[sol[2,1:300],sol[3,1:300]],label=["Exposed","Infected"],xlabel="Time [days]")
p2 = plot(sol.t,[sol[2,:],sol[3,:]],label=["Exposed","Infected"],xlabel="Time [days]")
plot(p1,p2, layout = l)
savefig("Wuhan_epidemic_curve_prediction.svg")

# Equilibrium point
tspan=(0.0,365.0*13)
I0v = [ 5 10 15 15 15 15 15 10 6 3 0 0 0 0]*1e4;
E0v = [ 0 0 0 1 2 3 4 4 4 4 3.5 2.5 1.5 0.5 ]*1e4;

l = @layout [a b]
p1 = plot(xlabel="Exposed",ylabel="Infected",legend=false)
p2 = plot(xlabel="Time [days]",ylabel="Exposed & Infected",xaxis=(:log10,(1,1e4)),legend=false)
for iv in 1:length(I0v)
   initial = [S, E0v[iv], I0v[iv], R, V]
   seirv_prob = ODEProblem(seirv_ode,initial,tspan,params)
   sol = solve(seirv_prob,alg=BS3(),saveat=1.0)
   plot!(p1,sol[2,:],sol[3,:])
   plot!(p2,sol.t,[sol[2,:],sol[3,:]])
end
plot(p1,p2,layout=l)
savefig("Wuhan_endemic_equilibrium.svg")


# Troubleshooting
# The right subfigure of Wuhan_endemic_equilibrium is blank when saved
# but if ploted interactively appears.
