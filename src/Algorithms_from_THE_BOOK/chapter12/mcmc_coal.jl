# Author: Kenneth Lange @ University of California, Los Angeles 

using StatsBase

"""Implements Gibb's sampling for a Poisson change-point model."""
function mcmc_coal(x::Vector{Int}, alpha::T, beta::T, 
  gamma::T, delta::T, trials::Int) where T <: Real
  n = length(x)
  lambda = mean(x) * ones(trials)
  mu = mean(x) * ones(trials)
  k = ones(Int, trials)
  kmass = zeros(n)
  for i = 2:trials
    ki = k[i - 1]
    lambda[i] = gamma_deviate(alpha + sum(x[1:ki]), beta + ki, 1)[1]
    mu[i] = gamma_deviate(gamma + sum(x[ki + 1:end]), delta + n - ki, 1)[1]
    for j = 1:n
      kmass[j] = exp(j * (mu[i] - lambda[i])) * (lambda[i] / mu[i])^sum(x[1:j])
    end
    k[i] = discrete_deviate(kmass, 1)[1]
  end
  lambda, mu, k
end

"""Generates n gamma deviates with shape parameter alpha
and intensity lambda."""
function gamma_deviate(alpha::T, lambda::T, n::Int) where T <: Real
#
  x = zeros(T, n)
  m = floor(Int, alpha)
  beta = alpha - m
  y, z = zero(T), zero(T)
  for i = 1:n
    z = - log(prod(rand(T, m)))
    if beta <= one(T) / 10^6
      y = zero(T)
    elseif beta < one(T) / 10^3
      y = (beta / rand(T))^(one(T) / (one(T) - beta))
    else
      r, s = one(T) / beta, beta - one(T)
      for trial = 1:1000
        u = rand(T, 2)
        y = - log(one(T) - u[1]^r)
        if u[2] <= (y / (one(T) - exp(- y)))^s
          exit
        end
      end
    end
    x[i] = (z + y) / lambda
  end
  x
end

"""Generates n random deviates according a given mass distribution."""
function discrete_deviate(mass::Vector{T}, n::Int) where T <: Real
  categories = length(mass)
  prob = mass / sum(mass)
  x = zeros(T, n)
  for trial = 1:n
    u = rand(T)
    s = zero(T)
    for i = 1:categories
      s = s + prob[i]
      if u <= s
        x[trial] = i
        break
      end
    end
  end
  x
end

alpha, beta, gamma, delta, trials = 4.0, 1.0, 1.0, 2.0, 10000
x = [4, 5, 4, 0, 1, 4, 3, 4, 0, 6, 3, 3, 4, 0, 2, 6,
  3, 3, 5, 4, 5, 3, 1, 4, 4, 1, 5, 5, 3, 4, 2, 5,
  2, 2, 3, 4, 2, 1, 3, 2, 2, 1, 1, 1, 1, 3, 0, 0,
  1, 0, 1, 1, 0, 0, 3, 1, 0, 3, 2, 2, 0, 1, 1, 1,
  0, 1, 0, 1, 0, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 2,
  3, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 2, 0, 0, 1, 4,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1];
lambda, mu, k = mcmc_coal(x, alpha, beta, gamma, delta, trials)
describe(k .+ 1850)
