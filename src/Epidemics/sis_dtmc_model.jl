
# SIS Discrete Time Markov Chain (DTMC) epidemic model
# =====================================================
# I(t): one independent random variable; S(t) = N - I(t)
# pᵢ(t) = Prob{I(t) = i}, ∀ i = 0,..., N and t = 0, Δt, 2Δt, ...; ∑ᵢ pᵢ(t) = 1
# Prob{I(t+Δt)|I(0),I(Δt),...,I(t)} = Prob{I(t+Δt)|I(t)} ← Markov property
# pⱼᵢ(t+Δt,t) = Prob{I(t+Δt) = j|I(t) = i} ← transition probability
# Here the deterministic model is autonomous, so the process is time homogeneous
# Therefore, pⱼᵢ(t+Δt,t) does not depend on time → pⱼᵢ(Δt).
# Δt sufficiently small such that the number of infected individuals changes by at most one.
#            /  βi(N-i)Δt/N,                  j = i + 1
#           |   (b + γ)iΔt,                   j = i - 1
# pⱼᵢ(Δt) = |   1 - [βi(N-i)/N + (b + γ)i]Δt, j = i
#            \  0,                            j ≠ i+1, i, i-1
# b(i) = βi(N-i)/N
# d(i) = (b + γ)i
# Δt chosen such that max{[b(i) + d(i)]Δt} ≤ 1; i ∈ {1,...,N}.
# pᵢ(t+Δt) = pᵢ₋₁(t)bt(i-1)Δt + pᵢ₊₁(t)d(i+1)Δt + pᵢ(t)(1 - [b(i) + d(i)]Δt).
# P(Δt) is the transition matrix with pⱼᵢ(Δt) elements; it is a stochastic matrix
# p(t+Δt) = P(Δt)p(t) = Pⁿ⁺¹(Δt)p(0) where t = nΔt


using Plots

Time = 2000
dtt = 0.01 # time step
β = 1dtt
b = 0.25dtt
γ = 0.25dtt
N = 100 # Total population size
en = 50 # plot every enth time interval
T = zeros(N+1,N+1) # T is the transition matrix, defined below
v = range(0,N,length=N+1)
p = zeros(Time+1,N+1)
p[1,3] = 1 # Two individuals initially infected: Prob{I(0)=2} = 1; R₀ = 2
bt = β*v .* (N .- v)/N
dt = (b + γ) .* v

for i in 2:N
   T[i,i]   = 1 - bt[i] - dt[i] # diagonal entries
   T[i,i+1] = dt[i+1]           # superdiagonal entries
   T[i+1,i] = bt[i]             # subdiagonal entries
end
T[1,1] = 1
T[1,2] = dt[2]
T[N+1,N+1] = 1 - dt[N+1]
for t in 1:Time
   y = T*p[t,:]
   p[t+1,:] = y'
end
#surface(p)
pm = zeros(div(Time,en)+1,N+1)
pm[1,:] = p[1,:]
for t in 1:div(Time,en)
   pm[t+1,:] = p[en*t,:]
end
#ti = range(0,Time,length=div(Time,en)+1)
#st = range(0,N,length=N+1)
wireframe(pm,xlabel="Infectives",ylabel="Time steps",zlabel="Probability")



