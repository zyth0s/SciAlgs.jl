
# https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/
# Stable periodic (Arenstorf) orbits around the Earth and the Moon
# were used to design a trajectory for the Apollo program

using DifferentialEquations

"""
Planar, circular restricted 3 body problem
xᵢ₊₁' = xᵢ' 
yᵢ₊₁' = yᵢ' 
xᵢ₊₁'' = xᵢ'  + 2yᵢ' - μ' (xᵢ+μ)/D1 - μ (xᵢ+μ')/D2
yᵢ₊₁'' = yᵢ'  - 2xᵢ' - μ' yᵢ    /D1 - μ yᵢ    /D2
with
D₁ = ((xᵢ+μ)² + yᵢ²)^(3/2)
D₂ = ((xᵢ-μp)² + yᵢ²)^(3/2)
"""
function planar_circ_3BP!(dv,v,params,t)
   μ = 0.01227747; μp = 1.0 - μ # Moon mass fraction, Earth mass fraction
   x, y, dx, dy = v                             # xᵢ, yᵢ, xᵢ', yᵢ'
   D1 = ((x+μ )^2 + y^2)^(3/2)
   D2 = ((x-μp)^2 + y^2)^(3/2)
   dv[1] = dx                                   # xᵢ₊₁'
   dv[2] = dy                                   # yᵢ₊₁'
   dv[3] = x + 2dy - μp*(x+μ)/D1 - μ*(x-μp)/D2  # xᵢ₊₁''
   dv[4] = y - 2dx - μp*y    /D1 - μ*y     /D2  # yᵢ₊₁''
end

# Arenstorf orbit with certain initial values
v0 = [0.994,0.0,0.0,-2.00158510637908] # x₀, y₀, x₀', y₀'
tspan = (0.0,17.06521656015796)

# Other Arenstorf orbit with certain initial values
#v0 = [1.2,0.0,0.0,-1.049357510] # x₀, y₀, x₀', y₀'

prob = ODEProblem(planar_circ_3BP!,v0,tspan,0)
v = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

using PyPlot

plot(v[1,:],v[2,:])
scatter(0,0)
scatter(1,0)
annotate(xy=(0.1,0.1),"Earth")
annotate(xy=(1.0,0.1),"Moon")

