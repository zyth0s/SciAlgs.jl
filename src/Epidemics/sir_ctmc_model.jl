
# http://epirecip.es/epicookbook/chapters/sir-stochastic-discretestate-continuoustime/julia
# Continuous Time Markov Chain (CTMC)
using DataFrames
using Distributions
import Random

function sir(beta,gamma,N,S0,I0,R0,tf)
    t = 0
    S = S0
    I = I0
    R = R0
    ta=Array{Union{Float64, Missing}}(undef,0)
    Sa=Array{Union{Float64, Missing}}(undef,0)
    Ia=Array{Union{Float64, Missing}}(undef,0)
    Ra=Array{Union{Float64, Missing}}(undef,0)
    while t < tf
        push!(ta,t)
        push!(Sa,S)
        push!(Ia,I)
        push!(Ra,R)
        pf1 = beta*S*I
        pf2 = gamma*I
        pf = pf1+pf2
        dt = rand(Exponential(1/pf))
        t = t+dt
        if t>tf
            break
        end
        ru = rand()
        if ru<(pf1/pf)
            S=S-1
            I=I+1
        else
            I=I-1
            R=R+1
        end
    end
    results = DataFrame()
    results[:time] = ta
    results[:S] = Sa
    results[:I] = Ia
    results[:R] = Ra
    return(results)
end


Random.seed!(42)

sir_out = sir(0.1/1000,0.05,1000,999,1,0,365);

head(sir_out)

using StatPlots

@df sir_out plot(:time, [:S :I :R], xlabel="Time",ylabel="Number")
