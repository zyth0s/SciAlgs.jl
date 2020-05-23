
# SIR model with age segregation
# Inspired by Sherry Towers R script
# Further info at http://sherrytowers.com/2012/12/11/sir-model-with-age-classes/

using DifferentialEquations
using LinearAlgebra: eigen
using Plots

function sir_age_ode!(dx,x,p,t)
   β, γ, C = p
   S = x[:,1]
   I = x[:,2]
   R = x[:,3]
   N = S+I+R
   dx[:,1] = -β*S .* (C * (I./N))        # dS
   dx[:,2] =  β*S .* (C * (I./N)) .- γ*I # dI
   dx[:,3] =  γ*I                        # dR
end

lcalculate_transmission_probability = 0 # if this is 1, then calculate the transmission probability from R0
                                        # otherwise, assume it is β=0.05

npop = 1e7
f = [0.25, 0.75] # two age groups, with 25% kids and 75% adults
N = npop*f
#--
nage = length(f)
I₀ = ones(nage)
S₀ = N .- I₀
R₀ = zeros(nage)
γ = 1/3.0 # average recovery time of 3 days
R0 = 1.5

# Contact matrix 
# number of contacts per day a person of a row group makes with members of a col group
#     kids   adults
C = [   18        9  # kids   ∴ have 18+9=27 contacts per day
         3       12] # adults ∴ have 3+12=15 contacts per day
# Remark: all kids have an adult in the home
# Remark: not all adults have kids
f[1]*C[1,2] == f[2]*C[2,1] || error("Contact matrix does not satisfy reciprocity")

β = 0.05
if lcalculate_transmission_probability == 1
   M = C
   for i in nage, j in nage
      M[i,j] = C[i,j] * f[i]/f[j]
   end
   eig = eigen(M)
   # reverse engineer β from the R0 and gamma 
   β = R0*γ/maximum(real(eig.values))  
end

parms = [β, γ, C]
initial = hcat(S₀, I₀, R₀)
tspan = (0.0,150.0)

sir_age_prob = ODEProblem(sir_age_ode!,initial,tspan,parms)
sir_age_sol = solve(sir_age_prob,saveat = 1);


println("The fraction of kids   that were infected is $(maximum(sir_age_sol[5,:])/N[1])")
println("The fraction of adults that were infected is $(maximum(sir_age_sol[6,:])/N[2])")
println("The total final size is $(maximum(sir_age_sol[5,:]+sir_age_sol[6,:])/npop)")

#plot(sir_age_sol,xlabel="Time",ylabel="Number", label=["Skids","Sadults","Ikids","Iadults","Rkids","Radults"])
plot(sir_age_sol.t, [sir_age_sol[3,:]/N[1],sir_age_sol[4,:]/N[2]],
     xlabel="Time, in days",
     ylabel="Fraction of each sub-population infected (prevalence)",
     title="Pandemic influenza simulation with age-structured SIR model",
     label=["kids","adults"])

