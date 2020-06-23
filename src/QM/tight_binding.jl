using LinearAlgebra
using Plots, Measures
using Parameters
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib
pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("sci")
mpl.use(backend="Qt5Agg")

include("dos_broadening.jl")

function buildH2(Δ,t)
  # He⁺
  [ 0 t; t Δ]
end
function buildHchain3(Δ1,Δ2,t12,t23)
  [  0 t12   0;
   t12  Δ1 t23;
     0 t23  Δ2]
end
function buildHring3(Δ1,Δ2,t12,t13,t23)
  [  0 t12 t13;
   t12  Δ1 t23;
   t13 t23  Δ2]
end
function buildHchain_tridiag(t,N=1001)
  H = SymTridiagonal(zeros(N),t*ones(N-1))
end
function buildHchain(t,N=1001,Δ=missing,impuritysite=missing)
  # homoatomic chain 
  if ismissing(Δ)
     return buildHchain_tridiag(t,N)
  end
  H = SymTridiagonal(zeros(N),t*ones(N-1))
  H = convert(Matrix,H)
  H[impuritysite,impuritysite] = Δ
  H
end
function buildHring(t,N=1001)
  H = SymTridiagonal(zeros(N),t*ones(N-1))
  H = convert(Matrix,H)
  H[N,1] = H[1,N] = t
  H
end
function buildHring_impurity(t,Δ,N=1001)
  H = SymTridiagonal(zeros(N),t*ones(N-1))
  H = convert(Matrix,H)
  H[N,1] = H[1,N] = t
  H[N÷2,N÷2] = Δ
  H
end
function tight_binding(hamiltonianConstructor::Function,params,name)
   @unpack t, N, Δ, impuritysite = params
   if ismissing(Δ)
      println("$name with $N sites with t = $t")
   else
      println("$name with $N sites with t = $t, impurity at $impuritysite, Δ = $Δ")
   end
   H = hamiltonianConstructor(params...)
   es, vs = eigen(H)
   W = 2abs(t) # half-bandwith 1D
   emargin = 1.5
   println("Lowest energy state with E = $(es[1])", es[1]<-W ? " (surface state)" : "")
   band = scatter(es, xlabel="Site", ylabel="Energy",ms=2,leg=false,
                  ylims=(-W-emargin,W+emargin))

   e_dos, dos = dos_broadening(es)
   dos = plot(dos,e_dos, xlabel="DOS",leg=false,
              xlims=(0,1),ylims=(-W-emargin,W+emargin))

   e_pdos1, pdos1 = dos_broadening(es,vs[1,:].*vs[1,:])
   pdos = plot(pdos1,e_pdos1,label="Site 1", ylims=(-W-emargin,W+emargin))
   e_pdos500, pdos500 = dos_broadening(es,vs[N÷2,:].*vs[N÷2,:])
   pdos = plot!(pdos, pdos500, e_pdos500, label="Site $(N÷2)",
                xlims=(0,0.0013),
                xlabel="PDOS")

   l = @layout [ a b c]
   plot(band,dos,pdos, layout=l, size=(1200,400),
        left_margin=5mm,bottom_margin=5mm,
        right_margin=0mm,top_margin=0mm)
   if ismissing(impuritysite)
      savefig("../figures/$(name)_t=$(t)_N=$(N)_E0=0.pdf")
   else
      savefig("../figures/$(name)_t=$(t)_N=$(N)_E0=0_@$(impuritysite)=$Δ.pdf")
   end
end

function tight_binding_1D_2sites(a,Δ₁₂,t)
   # 1D chain/ring with 2 sites per unit cell, s orbitals
   b = 2a # unit cell with two sites
   ε1 = 0 # without loss of generality energy reference at 0
   ε2 = Δ₁₂

   kpath = range(-π/b,stop=π/b,length=100)
   Enk = zeros(2,length(kpath))

   for (ik,k) in enumerate(kpath)
      H_k = [     ε1        2t*cos(k*a);
              2t*cos(k*a)       ε2      ]
      e, vs = eigen(H_k)
      Enk[:,ik] = e
   end

   plot(kpath/π,[Enk[1,:], Enk[2,:]],
        label=["Ground state" "Excited state"],
        xlabel="k pi/a", ylabel="Energy",
        leg = :inside,
       )
   savefig("../figures/1d_2sites.pdf")
end

function tight_binding_1D_2sites_sp_orbs(a,Δ₁₂,t)
   # 1D chain/ring with 2 sites per unit cell, sp orbitals
   # Redefined ψpk = i 1/√N * ∑R₂ exp(ikR₂) ϕp(r-R₂) to have a real matrix
   b = 2a # unit cell with two sites
   εs = Δ₁₂
   εp = 0 # without loss of generality energy reference at 0

   kpath = range(-π/b,stop=π/b,length=100)
   Enk = zeros(2,length(kpath))

   for (ik,k) in enumerate(kpath)
      H_k = [      εs        -2t*sin(k*a);
              -2t*sin(k*a)        εp      ]
      e, vs = eigen(H_k)
      Enk[:,ik] = e
   end

   plot(kpath/π,[Enk[1,:], Enk[2,:]],
        label=["p state" "s state"],
        xlabel="k pi/a", ylabel="Energy",
        leg = :right,
       )
   savefig("../figures/1d_2sites_sp_orbs.pdf")
end

function tight_binding_2D(a,b,Δ₁₂,ta,tb)
   # 2D rectangular lattice with 1 site per unit cell, s orbitals
   @assert a == b # does not handle rectangular with analytic formula
   @assert Δ₁₂ == 0 # does not handle ≠ ε with analytic formula
   ε1 = 0 # without loss of generality energy reference at 0
   ε2 = Δ₁₂

   Hamiltonian(kx,ky) = [               ε1               2ta*cos(kx*a) + 2tb*cos(ky*b);
                         2ta*cos(kx*a) + 2tb*cos(ky*b)                  ε2            ]

   # ............................
   # Band plot
   nkpts = 100
   kpath = zeros(2,3nkpts)
   kpath[1,  1:nkpts]  = range(0,stop=π/a,length=nkpts)
   kpath[2,  1:nkpts] .= 0
   #
   kpath[1,(nkpts+1):2nkpts] .= π/a
   kpath[2,(nkpts+1):2nkpts]  = range(0,stop=π/b,length=nkpts)
   #
   kpath[1,(2nkpts+1):3nkpts]  = range(π/a,stop=0,length=nkpts)
   kpath[2,(2nkpts+1):3nkpts]  = range(π/b,stop=0,length=nkpts)
   Enk = zeros(2,3nkpts)

   for ik in 1:3nkpts
      kx = kpath[1,ik]
      ky = kpath[2,ik]
      H_k = Hamiltonian(kx,ky)
      e, vs = eigen(H_k)
      #Enk[:,ik] = e # FIXME
      Enk[1,ik] = ε1 + ε2 - 2abs(ta)*cos(kx*a) - 2abs(tb)*cos(ky*b) # analytic formula
   end

   #plot(1:3nkpts,[Enk[1,:],Enk[2,:]],
   #plot(1:3nkpts,Enk[1,:],
   #     label=["Ground state" "Excited state"],
   #     xlabel="kpath/pi", ylabel="Energy",
   #     leg = :inside,
   #    )
   #savefig("2d_1sites.pdf")
   fig = pyplt.figure()
   fig.set_dpi(260)
   ax = fig.add_subplot(111)
   ax.plot(1:3nkpts,Enk[1,:], label="Ground state")
   ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(100))
   ax.set_xticklabels(["d","(0,0)T", "(1,0)T", "(1,1)T", "(0,0)T"])
   ax.set_xlim(0,3nkpts)
   ax.set_ylim(-2abs(ta)-2abs(tb),2abs(ta)+2abs(tb))
   ax.grid(true)
   ax.set_xlabel("(kx,ky)T pi/a")
   ax.set_ylabel("Energy")
   fig.tight_layout(pad=0.1)
   pyplt.savefig("../figures/2d_1sites.pdf")

   # ..........................................................
   # BZ integration
   nmesh = 101
   Enk = zeros(2,nmesh,nmesh)
   kxrange = range(-π/a,stop=π/a,length=nmesh)
   kyrange = range(-π/b,stop=π/b,length=nmesh)
   for (ikx,kx) in enumerate(kxrange), 
       (iky,ky) in enumerate(kyrange)
      H_k = Hamiltonian(kx,ky)
      #e, vs = eigen(H_k) # FIXME
      ##_,_,eschur = schur(H_k)
      Enk[1,ikx,iky] = ε1 + ε2 - 2abs(ta)*cos(kx*a) - 2abs(tb)*cos(ky*b) # analytic result a = b; ...
      #if abs(kx) + abs(ky) < pi/a
      #   Enk[:,ikx,iky] = e
      #else
      #   Enk[:,ikx,iky] = reverse(e)
      #end
      #isapprox(e[1] , ε1 + ε2 + 2ta*cos(kx*a) + 2tb*cos(ky*b), atol=1e-1) || error("wrong eigenvalue $(e[1]) ≠ $(ε1 + ε2 + 2ta*cos(kx*a) + 2tb*cos(ky*b))")
      #@show kx, ky
      #println(e[1]," ", e[2], " ", ε1 + ε2 - 2abs(ta)*cos(kx*a) - 2abs(tb)*cos(ky*b))
      #e[1], e[2] = e[2], e[1]
   end

   fig = pyplt.figure()
   fig.set_dpi(260)
   ax = fig.add_subplot(111)
   ax.set_aspect("equal")
   cont = ax.contourf(kxrange./π,kyrange./π,Enk[1,:,:], levels=15)
   #cont = ax.plot_surface(kxrange./π,kyrange./π,Enk[1,:,:])
   #cont = ax.imshow(Enk[1,:,:],extent=(-1,1,-1,1))
   # Set axis labels
   ax.set_xlabel("kx pi/a")
   ax.set_ylabel("ky pi/b")
   pyplt.colorbar(cont,format="%+3.1f")
   fig.tight_layout(pad=0.1)
   pyplt.savefig("../figures/2d_1sites_contour.pdf")
end

# -----------------------------------------------------------------------------------
# Two atomic orbitals at sites 1 and 2
t = 0; Δ = 2
println("A.1.  t = $t, Δ = $Δ")
H = buildH2(Δ,t)
es, vs = eigen(H)
for i in 1:length(es)
  println("Eigenstate with E = $(es[i]) and vector $(vs[i,:])")
end

# -----------------------------------------------------------------------------------
t = 2; Δ = 0
println("A.2. t = $t, Δ = $Δ")
H = buildH2(Δ,t)
es, vs = eigen(H)
for i in 1:length(es)
  println("Eigenstate with E = $(es[i]) and vector $(vs[i,:])")
end

# -----------------------------------------------------------------------------------
t = 1; Δ = 10
println("A.3.  t << Δ; t = $t, Δ = $Δ")
H = buildH2(Δ,t)
es, vs = eigen(H)
for i in 1:length(es)
  println("Eigenstate with E = $(es[i]) and vector $(vs[i,:])")
end

# -----------------------------------------------------------------------------------
# 3 centers 1 electron
println("B.1. Δ₁ = 2, Δ₂ = 5, t₁₂ = 1, t₁₃ = 2, t₂₃ = 1")
H = buildHring3(2,5,1,2,1)
es, vs = eigen(H)
for i in 1:length(es)
  println("Eigenstate with E = $(es[i]) and vector $(vs[i,:])")
end

# -----------------------------------------------------------------------------------
t = 2.5; N = 1001
println("C.1. Ring with $N sites and t = $t")
H = buildHring(t,N)
es, vs = eigen(H)
W = 2abs(t) # half-bandwith
band = scatter(es, xlabel="Site", ylabel="Energy",ms=2,leg=false,
               ylims=(-W-1,W+1))

e_dos, dos = dos_broadening(es)
dos = plot(dos,e_dos, ylabel="Energy", xlabel="DOS",leg=false,
            ylims=(-W-1,W+1))
l = @layout [ a b]
plot(band,dos,layout=l)
savefig("../figures/C.1_bandos_ring_t=$(t)_N=$(N)_E0=0.pdf")

# -----------------------------------------------------------------------------------
t = 2.5; N = 1001; Δ = -2; impuritysite = N÷2
println("C.2. Ring with $N sites with t = $t, impurity at $impuritysite Δ = $Δ")
H = buildHring_impurity(t,Δ,N)
es, vs = eigen(H)
println("Impurity state with E = $(es[1])")
W = 2abs(t) # half-bandwith
band = scatter(es, xlabel="Site", ylabel="Energy",ms=2,leg=false,
               ylims=(-W-1,W+1))
e_dos, dos = dos_broadening(es)
dos = plot(dos,e_dos, ylabel="Energy", xlabel="DOS",leg=false,
     ylims=(-W-1,W+1))
l = @layout [ a b]
plot(band,dos,layout=l)
savefig("../figures/C.2_bandos_ring_t=$(t)_N=$(N)_E0=0_middleimpurity=$Δ.pdf")

plot(vs[:,1].*vs[:,1],xlims=(490,510), label="Bound state",
     ylabel="Probability",xlabel="Site",
     title="Contribution of site orbitals to impurity state")
savefig("../figures/C.2_boundstate_ring_t=$(t)_N=$(N)_E0=0_middleimpurity=$Δ.pdf")


# Surface states: Schokley/Tamm states
## -----------------------------------------------------------------------------------
name = "D.1. Chain"
params = (t=2.5,N=1001,Δ=missing,impuritysite=missing)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
name = "D.2. Chain"
params = (t=2.5,N=1001,Δ=-2,impuritysite=500)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
name = "D.3. Chain"
params = (t=2.5,N=1001,Δ=-2,impuritysite=1)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
name = "D.4. Chain"
params = (t=2.5,N=1001,Δ=-4,impuritysite=1)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
name = "E.1. Finite Chain"
params = (t=2.5,N=9,Δ=-2,impuritysite=1)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
name = "E.2. Finite Chain"
params = (t=2.5,N=9,Δ=-4,impuritysite=1)
tight_binding(buildHchain,params,name)

## -----------------------------------------------------------------------------------
tight_binding_1D_2sites(1,2,2) # a, Δ₁₂, t

## -----------------------------------------------------------------------------------
tight_binding_1D_2sites_sp_orbs(1,2,2) # a, Δ₁₂, t

## -----------------------------------------------------------------------------------
tight_binding_2D(1,1,0,1,1) # a, b, Δ₁₂, ta, tb
#tight_binding_2D(10,5,2,0.5,1) # a, b, Δ₁₂, ta, tb
