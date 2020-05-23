
# Section 2.9 of Fundamentals of Astronomy by Karttunen et al.
# Perturbation of coordinates by refraction
# Calculation of the refraction angle, R,
# that tells how much higher the object seems to be compared to its
# true altitude, a.

include("units.jl")

using Plots

function refraction_correction(a,P,T)
   # a is the altitude in degrees
   # P is the pressure in hPa
   # T is the temperature in °C
   T = 273.15 + T # K
   if a > 15
      R = P/T * 0.00452tand(90-a)
   else
      R = P/T * (0.1594+0.0196a + 0.00002a^2)/(1 + 0.505a + 0.0845a^2)
   end
   R
end

plt = plot(title="Refraction of Earth's atmosphere with altitude",
           xlabel="Altitude, a [deg]",
           ylabel="Refraction angle, R [arcmin]")
for (P,T) in [(1050,-30),(950,30),(700,0)]
   Rarray = []
   aarray = range(0,stop=10,length=100) # in °
   for a in aarray
      R = refraction_correction(a,P,T)
      R = deg2min(R)
      push!(Rarray,R)
   end
   plot!(plt,aarray,Rarray,label="P = $P hPa, T = $T °C")
end
display(plt)
