
# Check gravity with a scale and weights.
# The weight of some object with mass m is
# modified as we add a mass madd to the Earth
# increasing its gravitational pull.

#   Initial              Final
#  
#    +---+               +---+
#    | m |     ^         | m |     ^    
#    +---+     |         +---+     |    
#    Scale     |         Scale     |
#      |       |           |       |    
#      |       |  F₀      madd     |  F₁
#      |       |           |       |    
#    +---+     |         +---+     |    
#    | Mₑ|     v         | Mₑ|     v    
#    +---+               +---+          
#  ASSUMPTION: the Earth radius at our position is the nominal Earth radius.
#  NOTE: we could use this to measure the Earth redius at our position.

using Printf

function gravitation_force_Newton(m₁,m₂,G,r₁₂)
   # G: gravitational constant
   # m₁, m₂: masses of the two objects
   # r₁₂: distance between objects
   G * m₁ / r₁₂^2 * m₂ 
end


function mass_from_force(F,a)
   # F: force
   # m: mass
   # a: acceleration
   # F = m * a² → m = F / a²
   F/a^2
end


G    = BigFloat(6.674e-11) # m³⋅kg⁻¹⋅s⁻²; gravitational constant
Mₑ   = BigFloat(5.9722e24) # kg; mass of the Earth (ₑ)
Rₑ   = BigFloat(6.3781e6)  # m; nominal Earth radius
m    = BigFloat(0.1)       # kg; mass of reference object
madd = BigFloat(10)        # kg; mass added behind the scale
scaleHeight = 1            # m; scale is placed 1 m above the floor

# Initial force
F₀ = gravitation_force_Newton(Mₑ,m,G,Rₑ+scaleHeight)

# Modified force
F₁ = gravitation_force_Newton(Mₑ+madd,m,G,Rₑ+scaleHeight)

ΔF = abs(F₁ - F₀)

Δm = mass_from_force(ΔF,G)
@printf("The scale shows a difference of %15.3f g \n",Δm*1e3)
@printf("when we add %15.3f kg to the Earth pull.",madd)
