
# Inspired by doi:10.1021/acs.jchemed.7b00003
#
using LinearAlgebra
#import REPL
using REPL.TerminalMenus
using Plots
unicodeplots()

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
   width,steps = p.width,p.steps
   xvec = range(-width/2,stop=width/2,length=steps)
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
   width, depth,steps = p.width,p.depth,p.steps
   xvec = range(-width,stop=width,length=steps)
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
   width1,width2,depth1,depth2,potwidth,steps = p.width1,p.width2,p.depth1,p.depth2,p.potwidth,p.steps
   A = 2(width1+width2+potwidth)
   xvec = range(-A,stop=A,length=steps)
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
   depth, ω, steps = p.depth, p.ω, p.steps
   width = sqrt(abs(2depth)/ω^2)
   A = 2width
   xvec = range(-A,stop=A,length=steps)
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
   ω,depth,steps = p.ω, p.depth,p.steps
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
   xvec = range(2start,stop=2stop,length=steps)
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
   depth,width,potwidth,nwells,steps = p.depth,p.width,p.potwidth,p.nwells,p.steps
   x_size = width * (div(nwells,2) + 0.5) + potwidth * (div(nwells,2) + 0.5)
   xvec = range(-x_size,stop=x_size,length=steps)
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

function schroedinger(p; periodic=false)
   xvec, U, Δx = potential(p)
   T = build_T(steps,Δx,periodic=periodic)
   E, V, nbound = diagonalize_hamiltonian(T,U)
   imprime(p,E,nbound)
   plot(xvec,U,
        xlabel = "x",
        ylabel = "Energy [Ha]",
        ylim = (minimum(U),maximum(U)),
        lab = "U",
        widen = false,
        #color = :red,
        width = 80,
        show=true)
   #V = V' * V # Ψ²
   plot(xvec,V[:,1],lab="Ground State Ψ₁")
   plot!(xvec,V[:,2], lab="Excited State Ψ₂", ylabel = "Prob.   ", widen=false, show=true)
   #xticks!(-8:8, ["$i" for i in -8:8])
   #xvec, U
end

function build_T(steps,Δx; periodic=false)
   ħ = 1.0
   m = 1.0
   Δ = SymTridiagonal(-2.0.*ones(steps),ones(steps-1))/Δx^2
   if periodic
      Δ = convert(Matrix,Δ)
      Δ[1,end] = Δ[end,1] = -1/Δx^2
   end
   -0.5*ħ^2/m * Δ
end

function diagonalize_hamiltonian(T,U)
   H = T + Diagonal(U)
   E,V = eigen(H)
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
if model == 1
   println("Enter the width of your infinite well in atomic units (a.u.) ∈ [0.5, 15]: ")
   width = parse(Float64,readline(stdin))
   println("Enter the number of wavefunctions you would like to plot ∈ ℕ: ")
   nwfn = parse(Int,readline(stdin))
   p = InfPotentialParams(width,nwfn,steps)
   schroedinger(p)
elseif model == 2
   println("Enter the width of your finite well in atomic units (a.u.) ∈ [1.0, 15]: ")
   width = parse(Float64,readline(stdin))
   println("Enter the depth of the well ∈ [20, 500]: ")
   depth = parse(Float64,readline(stdin))
   p = FinitePotentialParams(width, depth,steps)
   schroedinger(p)
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
   schroedinger(p)
elseif model == 4
   println("Enter the depth of the well ∈ [2, 15]: ")
   depth = parse(Float64,readline(stdin))
   println("Enter the force constant of your harmonic well ∈ [0.3, 1.4]. ")
   ω = parse(Float64,readline(stdin))
   p = HarmonicPotentialParams(depth, ω, steps)
   schroedinger(p)
elseif model == 5
   println("Enter the depth of the well ∈ [2, 15]: ")
   depth = parse(Float64,readline(stdin))
   println("Enter the force constant of your Morse well ∈ [0.05, 1.4]. ")
   ω = parse(Float64,readline(stdin))
   p = MorsePotentialParams(depth, ω, steps)
   schroedinger(p)
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
   schroedinger(p,periodic=true)
end

