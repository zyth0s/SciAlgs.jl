
# 1D Schroedinger solver
# Inspired by doi:10.1021/acs.jchemed.7b00003
#
using Parameters
using LinearAlgebra
#using Arpack
#using SparseArrays
#import REPL
using REPL.TerminalMenus
#using Plots
#unicodeplots(); #gr()
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib
pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("7-25in_color")
mpl.use(backend="Qt5Agg")

options = [
       "Particle in an infinite potential well",
       "Particle in a finite well",
       "Particle in a double finite well (unequal depth)",
       "Particle in a harmonic well",
       "Particle in a Morse well",
       "Kronig-Penney finite well"
   ]

menu = RadioMenu(options, pagesize=7)
model = request("Please enter the model potential you would like to study:",menu)

abstract type PESParams end

#---------------------------------------------
struct InfPotentialParams <: PESParams
   width; nwfn; steps
end

function imprime(params::InfPotentialParams,E,nbound)
   println("Well width (a.u.) $(params.width)")
   println("Lowest $(params.nwfn) states:")
   for i in 1:params.nwfn
      println("Energy: $(E[i])")
   end
end

function potential(p::InfPotentialParams)
   @unpack width,steps = p
   xvec = range(-width/2,width/2,length=steps)
   U = zeros(steps)
   Δx = xvec[2]-xvec[1]
   xvec, U, Δx
end
#---------------------------------------------
struct FinitePotentialParams <: PESParams
   width; depth; steps
end

function imprime(p::FinitePotentialParams,E,nbound)
   println("Well width (a.u.) $(p.width)")
   println("Well depth (a.u.) $(-p.depth)")
   println("Lowest $nbound bound states:")
   for i in 1:nbound
      println("Energy: $(E[i])")
   end
end

function potential(p::FinitePotentialParams)
   @unpack width, depth,steps = p
   xvec = range(-width,width,length=steps)
   Δx = xvec[2]-xvec[1]
   U = @. -depth*rect(xvec,width)
   xvec, U, Δx
end
#---------------------------------------------
struct DoubleSquarePotentialParams <: PESParams
   width1; width2
   depth1; depth2
   potwidth
   steps
end

function imprime(p::DoubleSquarePotentialParams,E,nbound)
   println("Well 1 width (a.u.) $(p.width1)")
   println("Well 2 width (a.u.) $(p.width2)")
   println("Well 1 depth (a.u.) $(-p.depth1)")
   println("Well 2 depth (a.u.) $(-p.depth2)")
   println("Potential width (a.u.) $(p.potwidth)")
   println("Lowest $nbound bound states:")
   for i in 1:nbound
      println("Energy: $(E[i])")
   end
end

function potential(p::DoubleSquarePotentialParams)
   @unpack width1,width2,depth1,depth2,potwidth,steps = p
   A = 2(width1+width2+potwidth)
   xvec = range(-A,A,length=steps)
   potwidth /= 2
   Δx = xvec[2]-xvec[1]
   U = @. -depth1*rect(xvec+potwidth+(width1/2),width1) -
           depth2*rect(xvec-potwidth-(width2/2),width2) 
   xvec, U, Δx
end
#---------------------------------------------
struct HarmonicPotentialParams <: PESParams
   depth; ω; steps
end

function imprime(p::HarmonicPotentialParams,E,nbound)
   println("Force constant $(p.ω)")
   println("Well depth (a.u.) $(abs(p.depth))")
   println("Lowest $nbound bound states:")
   for i in 1:nbound
      println("Energy: $(E[i])")
   end
end

function potential(p::HarmonicPotentialParams)
   @unpack depth, ω, steps = p
   width = sqrt(abs(2depth)/ω^2)
   A = 2width
   xvec = range(-A,A,length=steps)
   Δx = xvec[2]-xvec[1]
   U = @. 0.5ω^2 * xvec^2 - abs(depth)
   U[U .> 0] .= 0.0
   xvec, U, Δx
end
#---------------------------------------------
struct MorsePotentialParams <: PESParams
   depth; ω; steps
end

function imprime(p::MorsePotentialParams,E,nbound)
   println("Force constant $(p.ω)")
   println("Well depth (a.u.) $(abs(p.depth))")
   println("Lowest $nbound bound states:")
   for i in 1:nbound
      println("Energy: $(E[i])")
   end
end

function potential(p::MorsePotentialParams)
   @unpack ω,depth,steps = p
   depth = abs(depth)
   a = √(ω*depth/2)
   start = 0.0
   stop = 0.0
   morse_function(a,depth,x) =  depth * (exp(-2a*x) - 2exp(-a*x))
   while morse_function(a,depth,start) < 0.5depth
      start -= 0.01
   end
   while morse_function(a,depth,stop) < -0.1
      stop += 0.01
   end
   xvec = range(2start,2stop,length=steps)
   Δx = xvec[2]-xvec[1]
   U = morse_function.(a,depth,xvec)
   U[U .> 0] .= 0.0
   xvec, U, Δx
end
#---------------------------------------------
struct KronigPenneyPotentialParams <: PESParams
   depth
   width
   potwidth
   nwells
   steps
end

function imprime(p::KronigPenneyPotentialParams,E,nbound)
   println("Well width (a.u.) $(p.width)")
   println("Well depth (a.u.) $(abs(p.depth))")
   println("Well separation (a.u.) $(p.potwidth)")
   println("Lowest $nbound bound states:")
   for i in 1:nbound
      println("Energy: $(E[i])")
   end
end

function potential(p::KronigPenneyPotentialParams)
   @unpack depth,width,potwidth,nwells,steps = p
   x_size = width * (div(nwells,2) + 0.5) + potwidth * (div(nwells,2) + 0.5)
   xvec = range(-x_size,x_size,length=steps)
   Δx = xvec[2]-xvec[1]
   #
   U = @. -rect(xvec,width)
   for n in 1:(div(nwells,2) + 1)
      U -= @. rect(xvec + n*width + n*potwidth,width)
      U -= @. rect(xvec - n*width - n*potwidth,width)
   end
   U .*= depth
   xvec, U, Δx
end
#---------------------------------------------

function schroedinger(p; periodic=false,method=:centraldiff)
   xvec, U, Δx = potential(p)
   T = build_T(steps,Δx,periodic=periodic,method=method)
   E, V, nbound = diagonalize_hamiltonian(T,U)
   imprime(p,E,nbound)
   #plot(xvec,U,
   #     xlabel = "x",
   #     ylabel = "Energy [Ha]",
   #     ylim = extrema(U),
   #     lab = "U",
   #     widen = false,
   #     #color = :red,
   #     width = 80,
   #     show=true)
   #V = V' * V # Ψ²
   #plot(xvec,V[:,1],lab="Ground State Ψ₁")
   #plot!(xvec,V[:,2], lab="Excited State Ψ₂", ylabel = "Prob.   ", widen=false, show=true)
   #xticks!(-8:8, ["$i" for i in -8:8])
   #
   fig = pyplt.figure()
   #fig.set_size_inches([15/2.54,10/2.54],forward=true)
   fig.set_dpi(260)
   ax = fig.add_subplot(111)
   ax2 = ax.twinx()
   # Plot data
   ax.plot(xvec,U, label = "U",color=:black)
   ax2.plot(xvec,V[:,1],label="Ground state",color=:blue)
   ax2.plot(xvec,V[:,2],label="Excited state",color=:blue,":")
   # Set axes limits
   ax.set_xlim(extrema(xvec)...)
   ax.set_ylim(extrema(U)...)
   ax2.set_xlim(extrema(xvec)...)
   ax2.set_ylim(extrema(V[:,1:2])...)
   # Set axis labels
   ax.set_xlabel("x")
   ax.set_ylabel("Energy [Ha]")
   ax2.set_ylabel("Eigenstate amplitude", labelpad=10,color=:blue)
   # Edit the tick parameters of the new x-axis
   ax2.yaxis.set_tick_params(which="major", size=10, width=2, direction="in",labelcolor=:blue)
   ax2.yaxis.set_tick_params(which="minor", size=7, width=2, direction="in",labelcolor=:blue)
   # Add ticks manually to wfn axis
   ax2.yaxis.set_major_locator(mpl.ticker.FixedLocator(range(minimum(V[:,1:2]), maximum(V[:,1:2]), length=5)))
   ax2.yaxis.set_minor_locator(mpl.ticker.FixedLocator(range(minimum(V[:,1:2]), maximum(V[:,1:2]), length=19)))
   # Add tick labels manually to wfn axis
   #ax2.set_yticklabels(["-1.0", "-0.5","0.0", "0.5", "1.0"])
   # Add legend to plot
   #ax.legend(bbox_to_anchor=(1.0, 0.9), loc=1, frameon=false)
   #ax2.legend(bbox_to_anchor=(1.0, 0.7), loc=1, frameon=false)
   ax.legend(bbox_to_anchor=(-0.1, 1.15), loc="upper center", ncol=3,frameon=false)
   ax2.legend(bbox_to_anchor=(0.6, 1.15), loc="upper center", ncol=3,frameon=false)
   fig.tight_layout(pad=0.1) #,top_margin=3)
   #fig.subplots_adjust(top=0.15,left=0.15,right=0.2,bottom=0.1)
   #figsize = fig.get_size_inches()
   pyplt.draw()
end

function build_T(steps,Δx; periodic=false,method::Symbol)
   # atomic units (ħ = 1; m = 1) ⟹   ħ²/m = 1
   Δ = SymTridiagonal(fill(-2.0,steps),ones(steps-1))/Δx^2
   if periodic
      Δ = convert(Matrix,Δ)
      Δ[1,end] = Δ[end,1] = -1/Δx^2
   end
   if method == :numerov # doi:10.1119/1.4748813
      B = SymTridiagonal(fill(10.0,steps),ones(steps-1)) / 12.0
      if periodic 
         B = convert(Matrix,B)
         B[1,end] = B[end,1] = 1/12
      end
      Δ = inv(B) * Δ
   end
   -0.5 * Δ
end

function diagonalize_hamiltonian(T,U)
   H = T + Diagonal(U)
   #E, V = eigen(H)
   E, V = eigen(Matrix(H)) # FIXME for Numerov
   indices = sortperm(E)
   E = E[indices]
   V = V[:,indices]
   #if issparse(H)
   #   E,V = eigs(H,nev=2)
   #else
   #   E,V = eigen(H)
   #end
   nbound = 0
   while E[1+nbound] < 0
      nbound += 1
   end
   E, V, nbound
end


# Heaviside function                  0_____
heaviside(x) = 0.5*(1+sign(x)) # _____|               _____
# Rectangular function                        # _____|  0  |______
rect(center,width) = heaviside(center+width/2)-heaviside(center-width/2)
# Boxcar function
boxcar(x,center,with,A) = A*rect(x+center,width) # least ambiguous
#boxcar(x,a,b,A)         = A*(heaviside(x-a) - heaviside(x-b)) # clearest
#boxcar(x,startleft,width,A) # a bit ambiguous, better center

steps = 2000
method = :centraldiff
method = :numerov
function selector()
   if model == 1
      println("Enter the width of your infinite well in atomic units (a.u.) ∈ [0.5, 15]: ")
      width = parse(Float64,readline(stdin))
      println("Enter the number of wavefunctions you would like to plot ∈ ℕ: ")
      nwfn = parse(Int,readline(stdin))
      p = InfPotentialParams(width,nwfn,steps)
      schroedinger(p,method=method)
   elseif model == 2
      println("Enter the width of your finite well in atomic units (a.u.) ∈ [1.0, 15]: ")
      width = parse(Float64,readline(stdin))
      println("Enter the depth of the well ∈ [20, 500]: ")
      depth = parse(Float64,readline(stdin))
      p = FinitePotentialParams(width, depth,steps)
      schroedinger(p,method=method)
   elseif model == 3
      println("Enter the width of your first finite well in atomic units (a.u.) ∈ [0.5 and 10.0]: ")
      width1 = parse(Float64,readline(stdin))
      println("Enter the width of your second finite well in atomic units (a.u.) ∈ [0.5 and 10.0]: ")
      width2 = parse(Float64,readline(stdin))
      println("Enter the depth of the first well ∈ [30, 500]: ")
      depth1 = parse(Float64,readline(stdin))
      println("Enter the depth of the second well ∈ [30, 500]: ")
      depth2 = parse(Float64,readline(stdin))
      println("Enter the potential width ∈ [0.1, 10.0]: ")
      potwidth = parse(Float64,readline(stdin))
      p = DoubleSquarePotentialParams( width1,width2,depth1,depth2,potwidth,steps)
      schroedinger(p,method=method)
   elseif model == 4
      println("Enter the depth of the well ∈ [2, 15]: ")
      depth = parse(Float64,readline(stdin))
      println("Enter the force constant of your harmonic well ∈ [0.3, 1.4]. ")
      ω = parse(Float64,readline(stdin))
      p = HarmonicPotentialParams(depth, ω, steps)
      schroedinger(p,method=method)
   elseif model == 5
      println("Enter the depth of the well ∈ [2, 15]: ")
      depth = parse(Float64,readline(stdin))
      println("Enter the force constant of your Morse well ∈ [0.05, 1.4]. ")
      ω = parse(Float64,readline(stdin))
      p = MorsePotentialParams(depth, ω, steps)
      schroedinger(p,method=method)
   elseif model == 6
      println("Enter the width of wells in atomic units (a.u.) ∈ [1.0 and 15.0]: ")
      width = parse(Float64,readline(stdin))
      println("Enter the depth of the well ∈ [20, 500]: ")
      depth = parse(Float64,readline(stdin))
      println("Enter the potential width ∈ [1.0, 15.0]: ")
      potwidth = parse(Float64,readline(stdin))
      println("Enter the number of wells ∈ [3, 7]:")
      nwells = parse(Int,readline(stdin))
      p = KronigPenneyPotentialParams(depth,width,potwidth,nwells,steps)
      schroedinger(p,periodic=true,method=method)
   end
end

selector()
# TODO: Generalized Matrix Numerov method
