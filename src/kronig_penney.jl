# Kronig-Penney model
# ....................
# Julia implementation of the algorithm of
# https://github.com/shamim-hussain/kronig_penney.git
# by Shamim Hussain 

using Plots
#Constants
h  = 6.626e-34  # PlanckConstant, J s
c  = 2.998e8    # SpeedOfLightInVacuum, m / s
ħ  = 1.055e-34  # ReducedPlanckConstant, J s
mₑ = 9.109e-31  # ElectronMass, kg
e  = 1.602e-19  # ElementaryCharge, C

#User inputs
U_eV = 1       # eV
a    = 1.5e-10 # m
b    = 6e-10   # m

#Derived values
U = U_eV*e
α₀ = √(2*mₑ*U/(ħ^2))

f(ζ) = (1-2*ζ)/(2*√(ζ*(ζ-1)))*sin(a*α₀*√(ζ)) *sin(b*α₀*√(ζ-1)) +cos(a*α₀*√(ζ))*cos(b*α₀*√(ζ-1))

ζ = range(0.1, stop=10.0, length=100000)
ζ = Array{Complex,1}(ζ)
fζ = f.(ζ)
fζ[isnan.(fζ)] .= 1
ζ = Array{Float64,1}(ζ)
fζ = Array{Float64,1}(fζ)
#println(ζ)
#println(fζ)

# Plot cos k(a+b) vs zeta
plot(ζ,fζ,
     reuse=false,
     ylims=(minimum(fζ)-0.5,3),
     xlabel = "ζ (=E/U_0) →",
     ylabel = "f(ζ) (=RHS) →",
     title = "RHS of f(ζ) vs. ζ; U= $(U_eV) eV; a= $(a*1e10) Å and b= $(b*1e10) Å"
    )
plot!([ζ[1], ζ[end]], [1, 1],linecolor=:red)
plot!([ζ[1], ζ[end]], [-1, -1],linecolor=:green)
savefig("fzeta_vs_zeta.png")


function bandplots(ζ,fζ)
   flg = abs.(fζ) .<= 1
   pred = plot( xlabel = "Crystal momentum, k(radian/meter) →",
                ylabel = "Energy, E (eV) →",
                title = "Reduced zone representation of the E-k relationship; U= $(U_eV), eV; a= $(a*1e10) Å and b= $(b*1e10) Å",
                xrotation = 45,
                xticks = ([-π, -π/2, 0, π/2, π], ["-pi/(a+b)","-pi/(2(a+b))","0","pi/(2(a+b))","pi/(a+b)"])
             )
   pext = plot( xlabel = "Crystal momentum, k(radian/meter) →",
                ylabel = "Energy, E (eV) →",
                title = "Extended zone representation of the E-k relationship; U= $(U_eV), eV; a= $(a*1e10) Å and b= $(b*1e10) Å",
                xrotation = 45,
                xticks = ([-6*π, -5*π, -4*π, -3*π, -2*π, -π, 0, π, 2*π, 3*π, 4*π, 5*π, 6*π], 
                        vcat(["$(i)pi/(a+b)" for i in -6:-1],"0",["$(i)pi/(a+b)" for i in 1:6]))
             )
   plst=1; k=1
   while (!isempty(flg)) && (k<6)
       pos=findall(flg)
       if isempty(pos)
           break
       end
       pfst=plst+pos[1]-1
       flg=flg[pos[1]:end]
       pos=findall(map(!,flg))
       if isempty(pos)
           break
       end
       plst=pfst+pos[1]-1
       flg=flg[pos[1]:end]
       
       kv=acos.(fζ[pfst:plst-1]) #/(a+b)
       ev=ζ[pfst:plst-1]*U_eV
       if mod(k,2) == 0
           plot!(pred,[-reverse(kv), kv], [reverse(ev),ev], linecolor=:blue)
           plot!(pext,kv.-k*π, ev, linecolor=:red);
           plot!(pext,-reverse(kv).+k*π, reverse(ev), linecolor=:green);
       else
           plot!(pred,[kv, -reverse(kv)], [ev,reverse(ev)], linecolor=:blue)
           plot!(pext,kv.+(k-1)*π, ev, linecolor=:blue);
           plot!(pext,-reverse(kv).-(k-1)*π, reverse(ev), linecolor=:black)
       end
       k += 1
   end
   l = @layout [a{0.1428w} b]
   plot(pred,pext,layout=l,grid=true)
   savefig("kronig_penney_bands.png")
   #fig
end


bandplots(ζ,fζ)


#function reduced_bandplot(ζ,fζ)
#   flg = abs.(fζ) .<= 1
#   fig = plot(reuse=false,
#              xlabel = "Crystal momentum, k(radian/meter) →",
#              ylabel = "Energy, E (eV) →",
#              title = "Reduced zone representation of the E-k relationship for U= $(U_eV), eV; a= $(a*1e10) Å and b= $(b*1e10) Å",
#              xticks = ([-π, -π/2, 0, π/2, π], ["-π/(a+b)","-π/(2(a+b))","0","π/(2(a+b))","π/(a+b)"])
#             )
#   #grid on
#   plst=1; k=1
#   while (!isempty(flg)) && (k<6)
#       pos=findall(flg)
#       if isempty(pos)
#           break
#       end
#       pfst=plst+pos[1]-1
#       flg=flg[pos[1]:end]
#       pos=findall(map(!,flg))
#       if isempty(pos)
#           break
#       end
#       plst=pfst+pos[1]-1
#       flg=flg[pos[1]:end]
#       
#       kv=acos.(fζ[pfst:plst-1]) #/(a+b)
#       ev=ζ[pfst:plst-1]*U_eV
#       if mod(k,2) == 0
#           plot!([-reverse(kv), kv], [reverse(ev),ev], linecolor=:blue)
#       else
#           plot!( [kv, -reverse(kv)], [ev,reverse(ev)], linecolor=:blue)
#       end
#       k += 1
#   end
#   savefig("reduced_band.png")
#   #fig
#end
#
#function extended_bandplot(ζ,fζ)
#   flg = abs.(fζ) .<= 1
#   fig = plot(reuse=false,
#              xlabel = "Crystal momentum, k(radian/meter) →",
#              ylabel = "Energy, E (eV) →",
#              title = "Extended zone representation of the E-k relationship for U= $(U_eV), eV; a= $(a*1e10) Å and b= $(b*1e10) Å",
#              xticks = ([-6*π, -5*π, -4*π, -3*π, -2*π, -π, 0, π, 2*π, 3*π, 4*π, 5*π, 6*π], 
#                        vcat(["$(i)π/(a+b)" for i in -6:-1],"0",["$(i)π/(a+b)" for i in 1:6]))
#             )
#   #xtickangle(45)
#   #grid on
#   plst=1; k=1
#   while (!isempty(flg)) && (k<6)
#       pos=findall(flg);
#       if isempty(pos)
#           break
#       end
#       pfst=plst+pos[1]-1;
#       flg=flg[pos[1]:end];
#       pos=findall(map(!,flg));
#       if isempty(pos)
#           break
#       end
#       plst=pfst+pos[1]-1;
#       flg=flg[pos[1]:end];
#       
#       kv=acos.(fζ[pfst:plst-1]) #/(a+b);
#       ev=ζ[pfst:plst-1]*U_eV;
#       if mod(k,2) == 0
#          plot!(kv.-k*π, ev, linecolor=:red);
#          plot!(-reverse(kv).+k*π, reverse(ev), linecolor=:green);
#       else
#          plot!(kv.+(k-1)*π, ev, linecolor=:blue);
#          plot!(-reverse(kv).-(k-1)*π, reverse(ev), linecolor=:black)
#       end
#       k += 1;
#   end
#   savefig("extended_band.png")
#   #fig
#end
#reduced_bandplot(ζ,fζ)
#extended_bandplot(ζ,fζ)
