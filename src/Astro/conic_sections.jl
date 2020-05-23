
# Conic sections in polar coordinates
# r = p/(1+ecosθ) with p = semi-latus rectum
# Circle:    e = 0     , p = a
# Ellipse:   e ∈ (0,1) , p = a(1-e²)
# Parabola:  e = 1     , p = 1/2a
# Hyperbola: e > 1     , p = a(e²-1)
# a is chosen to align all curves
# see Section A.2 in Fundamental Astronomy by Karttunen

using Plots

plt = plot(title="Conic sections",leg=false)
#for e in [0.0,0.6,1,1.4]
for e in range(0, stop=1.4, length=20)
   a = 1.0 # radius of the circle
   if e ≈ 0
      p = a*(1-e^2)
   elseif e < 1
      a = a/(1-e) # align with circle
      p = a*(1-e^2)
   elseif e ≈ 1
      a = 1/(4a) # align with circle
      p = 1/(2a)
   elseif e > 1
      a = a/(e-1) # align with circle
      p = a*(e^2-1)
   end
   Xarray = []
   Yarray = []
   for θ in range(0,stop = 2π, length=100)
      r = p/(1+e*cos(θ))
      push!(Xarray,r*cos(θ))
      push!(Yarray,r*sin(θ))
   end
   plot!(plt,Xarray,Yarray, label="e = $e")
end
xlims!(-4,4)
ylims!(-4,4)
display(plt)

