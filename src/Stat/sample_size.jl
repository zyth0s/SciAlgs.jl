
# Estimate the sample size needed to achieve a fixed
# maximum error when the population size is known

function Z(α)
   #Valor de Zα      	1.28	1.65	1.69	1.75	1.81	1.88	1.96
   #Nivel de confianza	80%	90%	91%	92%	93%	94%	95%
   if α == 0.05     # 95% CI
      return 1.96
   elseif α == 0.06 # 94% CI
      return 1.88
   elseif α == 0.07 # 93% CI
      return 1.81
   elseif α == 0.08 # 92% CI
      return 1.75
   elseif α == 0.09 # 91% CI
      return 1.69
   elseif α == 0.10 # 90% CI
      return 1.65
   elseif α == 0.20 # 80% CI
      return 1.28
   end
end

function sample_size(N,Zα, p, e)
    N * Zα^2 * p * (1-p) / ( e^2 * (N - 1) + Zα^2 * p * (1-p) )
end

N  = 47e6 # population size
e  = 0.01 # maximum error admited
α  = 0.05 # The confidence interval is 95% = (1-α) -> α = 0.05
p  = 0.5  # probability of the given characteristic; unknown but assumed; maximizes p(1-p)
Zα = Z(α) # normal distribution value that depends on α (confidence interval)

n = sample_size(N,Zα,p,e)
println("The required sample size is $n.")
