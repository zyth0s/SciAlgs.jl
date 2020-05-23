
# SIS Continuous Time Markov Chain (CTMC) epidemic model
# Three sample paths and the deterministic solution

using Random
using Plots

Random.seed!(42)

β = 1
b = 0.25
γ = 0.25
N = 100 # Total population size
init = 2
Time = 25
sim = 3
plot(xlabel="Time")
for k in 1:sim
   t=Array{Float64}(undef,0)
   s=Array{Float64}(undef,0)
   i=Array{Float64}(undef,0)
   push!(t,0)
   push!(i,init)
   push!(s,N-init)
   j = 1
   #while t[j] < Time
   while (i[j] > 0) & (t[j] < Time)
      u1 = rand() # uniform random number
      u2 = rand() # uniform random number
      a = (β/N) * i[j] * s[j] + (b + γ) * i[j]
      probi = (β*s[j]/N)/(β*s[j]/N+b+γ)
      push!(t,t[j] - log(u1)/a)
      if u2 <= probi
         push!(i, i[j]+1)
         push!(s, s[j]-1)
      else
         push!(i, i[j]-1)
         push!(s, s[j]+1)
      end
      j = j+1
   end
   plot!(t,i,label="Sample $k")
end

gui()
