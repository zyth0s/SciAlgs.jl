
# Notation
# a is the semimajor axis [AU]
# e is the eccentricity
# I is the inclination [arcdeg]
# Ω is the longitude of the ascending node [arcdeg]
# ϖ is the longitude of perihelion [arcdeg]
# λ is the mean longitude [arcdeg]
# M is the mean anomaly [arcdeg]
# E is the eccentric anomaly [arcdeg]
# f is the true anomaly [arcdeg or rad]
# ω is the argument of the perihelion [arcdeg]
# Newton's dot notation was used to denote variation of parameters with time.
# TODO: https://ssd.jpl.nasa.gov/?planet_pos

include("time.jl")
include("units.jl")

import Roots
using Plots
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib
pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("5in_color")
mpl.use(backend="Qt5Agg")

function orbital_elements_planet(t, T, object="EarthMoon")
   # t is Julian date since J2000.0 epoch
   # T is Julian century since J2000.0 epoch
   # object is the planet or planet+satellite
   # Accuracty ∼ 1 arcmin
   #  [a] = adim,[e] = adim,[I] = deg,[Ω] = deg,[ϖ] = deg,[λ] = deg
   # Taken from Table C.12 of Fundamental Astronomy 6ed
   if object == "Mercury"
      a = 0.38709893 + 0.00000066T
      e = 0.20563069 + 0.00002527T
      I = 7.00487 - sec2deg(23.51)*T
      Ω = 48.33167 - sec2deg(446.30)*T
      ϖ = 77.45645 + sec2deg(573.57)*T
      λ = 252.25084 + 4.09233880t
      return a,e,I,Ω,ϖ,λ
   elseif object == "Jupiter"
      a = 5.20336301 + 0.00060737T
      e = 0.04839266 - 0.00012880T
      I = 1.30530 - sec2deg(4.15)*T
      Ω = 100.55615 + sec2deg(1217.17)*T
      ϖ = 14.75385 + sec2deg(839.93)*T
      λ = 34.40438 + 0.08308676t
      return a,e,I,Ω,ϖ,λ
   end
end

# Example 6.1 Find the orbital elements of Jupiter on August 23, 1996.
t = julian_date(1996,8,23)
t_J2000 = days_since_J2000(t)
T = julian_century(t)
a,e,I,Ω,ϖ,λ = orbital_elements_planet(t_J2000,T,"Jupiter")
a ≈ 5.2033
e ≈ 0.0484
I ≈ 1.3053 # deg
Ω ≈ 100.5448 # deg
ϖ ≈ 14.7460 # deg
mod(λ,360) ≈ 292.540 # deg

#------------------------------------------------------------------

function planet_elements(planet,jdate)
   # Data taken from Solar System Dynamics Appendix A.4
   # a = a₀ + ȧt              [AU]
   # e = e₀ + ėt
   # I = I₀ + (İ/3600)t       [arcdeg]
   # ϖ = ϖ₀ + (ϖ̇/3600)t       [arcdeg]
   # Ω = Ω₀ + (Ω̇/3600)t       [arcdeg]
   # λ = λ₀ + (λ̇/3600+360Nᵣ)t [arcdeg]
   # where t is the time in Julian centuries since JD2000.0
   # Accuracy: max(error) = 600 arcsec in the period 1800-2050
   id = Dict("Mercury" => 1, "Venus" => 2, "Earth" => 3, "Mars" => 4,
             "Jupiter" => 5, "Saturn" => 6, "Uranus" => 7, "Neptune" => 8,
             "Pluto" => 9)[planet]
   a₀ = [0.38709893,     0.72333199,     1.00000011, 
                  1.52366231,     5.20336301,     9.53707032,
                  19.19126393,    30.06896348,    39.48168677]
   ȧ  = [0.00000066,     0.00000092,    -0.00000005,
                 -0.00007221,     0.00060737,    -0.00301530,
                 0.00152025,    -0.00125196,    -0.00076912]
   e₀ = [0.20563069,     0.00677323,     0.01671022,
                  0.09341233,     0.04839266,     0.05415060,
                  0.04716771,     0.00858587,     0.24880766]
   ė  = [0.00002527,    -0.00004938,    -0.00003804,
                  0.00011902,    -0.00012880,    -0.00036762,
                  -0.00019150,     0.00002514,     0.00006465]
   I₀ = [7.00487,        3.39471,        0.00005,
                  1.85061,        1.30530,        2.48446,
                  0.76986,        1.76917,       17.14175]
   İ  = [-23.51,          -2.86,         -46.94,
                -25.47,          -4.15,           6.11,
                -2.09,          -3.64,          11.07]
   ϖ₀ = [77.45645,      131.53298,      102.94719,
                336.04084,       14.75385,       92.43194,
                170.96424,       44.97135,      224.06676]
   ϖ̇  = [573.57,        -108.80,        1198.28,
               1560.78,         839.93,       -1948.89,
               1312.56,        -844.43,        -132.25]
   Ω₀ = [48.33167,       76.68069,      -11.26064,
                 49.57854,      100.55615,      113.71504,
                 74.22988,      131.72169,      110.30347]
   Ω̇  = [-446.30,        -996.89,      -18228.25,
              -1020.19,        1217.17,       -1591.05,
              1681.40,        -151.25,         -37.33]
   λ₀ = [252.25084,      181.97973,      100.46435,
                355.45332,       34.40438,       49.94432,
                313.23218,      304.88003,      238.92881]
   λ̇  = [261628.29,      712136.06,     1293740.63,
              217103.78,      557078.35,      513052.95,
              246547.79,      786449.21,      522747.90]
   Nᵣ = [415,            162,             99,
                 53,              8,              3,
                 1,              0,              0]
   t = julian_century(jdate)
   a = a₀[id] + ȧ[id]*t             
   e = e₀[id] + ė[id]*t
   I = I₀[id] + İ[id]*t/3600
   ϖ = ϖ₀[id] + ϖ̇[id]*t/3600
   Ω = Ω₀[id] + Ω̇[id]*t/3600
   λ = λ₀[id] + (λ̇[id]/3600+360Nᵣ[id])*t
   I = mod(I,360); ϖ = mod(ϖ,360); Ω = mod(Ω,360)
   λ = mod(λ,360)
   a,e,I,ϖ,Ω,λ
end

function orbit_xyz(planet,jdate,npts)
   a,e,I,ϖ,Ω,λ = planet_elements(planet,jdate)
   ω = mod(ϖ - Ω, 360)
   f = range(0,stop=2π,length=npts)
   cf = cos.(f)
   sf = sin.(f)
   sl = a*(1-e^2) # ellipse semi-latus
   r = map(f -> sl/(1+e*cos(f)),f) # polar ellipse
   x = r.*cf
   y = r.*sf
   z = zeros(npts)
   p = euler_rot(Ω,I,ω)
   mapslices( x -> p*x, vcat(x',y',z'),dims=[1])
end

function orbit_position(planet,jdate)
   a,e,I,ϖ,Ω,λ = planet_elements(planet,jdate)
   ω = mod(ϖ - Ω, 360)
   M = mod(λ - ϖ,360)
   #eq_Kepler(E) = E - e*sin(E) - deg2rad(M) # = 0
   #E = rad2deg(Roots.find_zero(eq_Kepler, deg2rad(M)))
   eq_Kepler(E) = E - rad2deg(e)*sind(E) - M # = 0
   E = Roots.find_zero(eq_Kepler, M)
   x = a*(cosd(E) - e)
   y = a*sqrt(1-e^2)*sind(E)
   p = euler_rot(Ω,I,ω)
   p*[x;y;0]
end

function euler_rot(Ω,I,ω)
   p1 = [ cosd(ω)    -sind(ω)          0;
          sind(ω)     cosd(ω)          0;
               0          0          1 ]
   p2 = [      1          0          0;
               0     cosd(I)    -sind(I);
               0     sind(I)     cosd(I) ]
   p3 = [ cosd(Ω)    -sind(Ω)          0;
          sind(Ω)     cosd(Ω)          0;
               0          0          1 ] 
   p1*p2*p3
end

function diagram_solar_system(Y::Int, m::Int, d::Int, ut)
   # Plot planets projected onto the ecliptic of J2000.0
   # at any datetime
   jd = julian_date_Montenbruck1989(Y,m,d,ut)
   default(markerstrokewidth=0,background_color=:black)
   planet_colors = Dict("Mercury" => :grey,"Venus" => :goldenrod, "Earth" => :mediumblue,"Mars" => :red,
                        "Jupiter" => :crimson,"Saturn" => :darkkhaki,"Uranus" => :mediumturquoise,"Neptune" => :royalblue,
                        "Pluto" => :darksalmon
                       )
   inner_planets = ["Mercury","Venus", "Earth","Mars"]
   outer_planets = ["Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]

   fig, axes = fig = pyplt.subplots(ncols=2,nrows=1)
   fig.set_dpi(200)
   #fig.suptitle("Solar system snapshot at $Y-$m-$d-$(trunc(Int,ut))h")
   axes[1].set_xlabel("Distance from Sun [AU]")
   axes[1].set_ylabel("Distance from Sun [AU]")
   axes[2].set_xlabel("Distance from Sun [AU]")
   #axes[2].set_ylabel("Distance from Sun [AU]")
   axes[1].set_title("Inner planets",y=1.04)
   axes[2].set_title("Outer planets",y=1.04)
   axes[1].scatter(0,0,color=:yellow, label="Sun", s=5)
   axes[2].scatter(0,0,color=:yellow, label="Sun", s=5)
   axes[1].set_xlim(-1.7,1.7)
   axes[1].set_ylim(-1.7,1.7)
   axes[2].set_xlim(-46,46)
   axes[2].set_ylim(-46,46)
   #for spine in ["left","top","right","bottom"]
   #   axes[1].spines[spine].set_visible(false)
   #   axes[2].spines[spine].set_visible(false)
   #end
   inner = scatter([0],[0],color=:yellow, label="Sun",
               markersize=5,
               xlabel="Distance from Sun [AU]",
               ylabel="Distance from Sun [AU]",
               leg = false,
               xlims = (-1.7,1.7), ylims=(-1.7,1.7),
               title="Inner planets")
   for planet in inner_planets
      a,e,I,ϖ,Ω,λ =planet_elements(planet,jd)
      xy = orbit_xyz(planet,jd,100)[1:2,:]
      plot!(inner,xy[1,:],xy[2,:],color=planet_colors[planet],label=nothing,
            line = (:solid,2))
      x,y,_ = orbit_position(planet,jd)
      scatter!(inner,[x],[y],label=nothing,color=planet_colors[planet],
               markersize=7)
      annotate!(inner,x+0.05,y, Plots.text(planet[1],planet_colors[planet],:left,15))
      axes[1].plot(xy[1,:],xy[2,:],color=planet_colors[planet],label=nothing,
                   linestyle=:solid,linewidth=2)
      axes[1].scatter(x,y,label=nothing,color=planet_colors[planet], s = 7)
      axes[1].annotate(string(planet[1]),xy=(x+0.05,y),color=planet_colors[planet])
   end

   outer = scatter([0],[0],color=:yellow, label="Sun",
               markersize = 5,
               xlabel = "Distance from Sun [AU]",
               #ylabel="Distance from Sun [AU]",
               leg = false,
               xlims = (-46,46), ylims=(-46,46),
               title = "Outer planets")
   for planet in outer_planets
      a,e,I,ϖ,Ω,λ =planet_elements(planet,jd)
      xy = orbit_xyz(planet,jd,100)[1:2,:]

      plot!(outer,xy[1,:],xy[2,:],color=planet_colors[planet],label=nothing,
            line = (:solid,2))
      x,y,_ = orbit_position(planet,jd)
      scatter!(outer,[x],[y],label=nothing,color=planet_colors[planet],
               markersize=7)
      annotate!(outer,x+1.5,y, Plots.text(planet[1],planet_colors[planet],:left,15))

      axes[2].plot(xy[1,:],xy[2,:],color=planet_colors[planet],label=nothing,
                   linestyle=:solid,linewidth=2)
      axes[2].scatter(x,y,label=nothing,color=planet_colors[planet], s = 7)
      axes[2].annotate(string(planet[1]),xy=(x+1.5,y),color=planet_colors[planet])
   end
   l = @layout [a b]
   plot(inner,outer,layout = l,size=(1200,650))
   savefig("solar_system_$Y-$m-$d-$ut.pdf")
   fig.tight_layout(pad=0.1) # top_margin conflicts with savefig
   pyplt.savefig("../figures/solar_system_$Y-$m-$d-$ut.pdf")
   pyplt.close("all")
end
# -----------------------------------------------------------
#1993 Sept 25 17ʰ32ᵐ  British Summer Time (GMT+1)
jd = julian_date_Montenbruck1989(1993,9,25,16+32.0/60.0)

planet_elements("Earth",jd)

# -----------------------------------------------------------
#1996 Aug 23 5ʰ32ᵐ British Summer Time (GMT+1)
diagram_solar_system(1993,9,25,16+32.0/60.0)


