
# http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-discretetime/julia
# This is a chain binomial model, a well-known Discrete Time Markov Chain (DTMC) model
# At each step:
# I(t+1) ∼ I(t) + Binom(S(t), ifrac)
# S(t+1) = S(t) - I(t+1)
# R(t+1) ∼ R(t) + Binom(I(t), rfrac)
using RandomNumbers
using DataFrames

# Binomial distribution with (n, p) parameters and rng seed
@inline @fastmath function randbn(n,p,rng)
    q = 1.0 - p
    s = p/q
    a = (n+1)*s
    r = exp(n*log(q))
    x = 0
    u = rand(rng)
    while true
        if (u < r)
            return x
        end
        u -= r
        x += 1
        r *= (a/x)-s
    end
end

# Step in the evolution
@inline @fastmath function sir(u, parms, rng)
    (S, I, R, Y) = u
    (β, γ, ι, N, δt) = parms
    λ = β * (I + ι) / N # Rate of contracting the disease per infective per δt
    ifrac = 1.0 - exp(-λ * δt) # Probability that the susceptible gets infected
    rfrac = 1.0 - exp(-γ * δt) # Probability that the infected recovers
    infection = randbn(S, ifrac, rng)
    recovery = randbn(I, rfrac, rng)
    return (S - infection, I + infection - recovery, R + recovery, Y + infection)
end

function simulate(r)
    parms = (0.1, 0.05, 0.01, 1000.0, 0.1)
    tf = 200
    t = 0:0.1:tf
    tl = length(t)
    S = zeros(tl)
    I = zeros(tl)
    R = zeros(tl)
    Y = zeros(tl)
    u0 = (999, 1, 0, 0)
    (S[1],I[1],R[1],Y[1]) = u0
    u = u0
    for i in 2:tl
        u = sir(u, parms, r)
        (S[i],I[i],R[i],Y[i]) = u
    end
    return DataFrame(Time=t,S=S,I=I,R=R,Y=Y)
end

seed = 42
r = Xorshifts.Xorshift128Plus(seed);

sir_out = simulate(r);

head(sir_out)

using Plots, StatPlots

@df sir_out plot(:Time, [:S :I :R], colour = [:red :green :blue], xlabel="Time",ylabel="Number")


