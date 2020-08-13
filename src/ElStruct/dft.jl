# A simple plane wave Density Functional Theory (DFT) code
# https://github.com/mpanho/simple-DFT/blob/master/DFT_1D.ipynb
import FFTW: fft, ifft
import LinearAlgebra: eigen, Diagonal, diagm #, ishermitian, LAPACK.syev!
import SpecialFunctions: erfcx
import Interpolations: CubicSplineInterpolation

function FFT(x)
  tmp = fft(x) / nG
  if maximum(abs.(imag.(tmp))) < 1e-10
    return real.(tmp)
  else
    return tmp
  end
end

function IFFT(x)
  tmp = ifft(x) * nG
  if maximum(abs.(imag.(tmp))) < 1e-10
    return real.(tmp)
  else
    return tmp
  end
end

function kinetic(k)
  # The kinetic operator is diagonal
  Diagonal((k.+ G).^2 ./rs^2)
end

function potential_op(V)
  tmp = zeros(nG,nG)
  #tmp += Diagonal(ones(nG)*V[1])
  tmp += diagm( 0 => ones(nG)*V[1])
  for i in 0:(nG-2)
    tmp += diagm( i+1 => ones(nG-i-1)*V[end-i])
    tmp += diagm(-i-1 => ones(nG-i-1)*V[i+2])
  end
  tmp
end

function erwin(V)
  bandst = zeros(nk,nB)
  wfk = zeros(Complex{Float64},nk,nG,nB)
  for i in 1:nk
    H = kinetic(kvec[i])+V
    #@assert ishermitian(H)
    #eig,v = eigen(H)
    eig, v = LAPACK.syev!('V','L',H)
    bandst[i,:] = eig[1:nB]
    wfk[i,:,:] = v[:,1:nB]
  end
  bandst, wfk
end

function occupation(eigenvalue)
  if (eigenvalue-EFermi)/kT > 20
    return 0.0
  elseif (eigenvalue-EFermi)/kT < -20
    return 1.0
  else #if -20 < (eigenvalue - EFermi) < 20
    return 1/(exp((eigenvalue-EFermi)/kT)+1)
  end
end

function calc_ρ(occ,wfk)
  ρ = zeros(nG)
  for ik in 1:nk, jn in 1:nB
    ρ += 2occ[ik,jn] .* abs.(IFFT(wfk[ik,:,jn])).^2
  end
  ρ /= 2a*nk
end

function coul_pot(r,rs)
  # 1D wire, harmonic trap
  # erfcx(x) = e^(x^2) erfc(x)
  # V(x) = √π/b * e^(x^2/4b^4) * erfc(|x|/2b)
  # eq. 4 in Ground state properties of the one dimensional Coulomb gas.
  # Casula, Sorella, ...
  b = 0.1/rs
  @. erfcx(abs(r)/(2b)) * √π/(b*rs)
end

function coul_per(r)
  coul_pot(r,rs) + coul_pot(r+a*nk,rs) + coul_pot(r-a*nk,rs)
end

function Ecorr(rs)
  # ≈ eq. 4 in "Density Functional theory beyond the linear regime: 
  # Validating adiabatic LDA"
  A = 4.66
  B = 2.092
  C = 3.735
  n = 1.379
  α = 23.63
  β = 109.9
  m = 1.837
  -(rs/(A + B*rs^n + C*rs^2)) * log(1 + α*rs + β*rs^m)
end

function rs_ρ(ρ)
  rs/(2ρ)
end

function VH(ρ)
  res = zeros(nG)
  for ix in 1:nG, ixp in -nG*nk/2:nG*nk/2
    res[ix] += ρ[Int(mod(ix-ixp,nG))+1] * coul_per(ixp*a/nG)
  end
  res * a / nG
end

# ----------------------------------------------------
# Set the parameters
Ne = 1 # number of electrons in the unit cell
rs = 1 # average r_s value
nG = 20 # number of plane-wave G vectors
nk = 40 # number of k vectors
@assert mod(nk,4) == 0
nB = 10 # number of bands
kT = 0.001 # Smearing parameter for a Fermi-Dirac distribution in Ry units

a = 2Ne # length of the unit cell
Nes = Ne*nk # number of electrons in the super-cell
kF = π/4 # Fermi wave vector corresponding to the mean density of the system
ρ0 = Ne/a # mean density
ΔG = 2π/a # reciprocal lattice spacing
x = range(0,a-1e-7,step=a/nG) # real space grid
kvec = range(-ΔG*(1/2-1/nk),ΔG/2+0.000001,length=nk) # BZ1 discretization ∈ (-2π/a,2π/a)
Gtmp = 0:1:(nG-1)
G = @. (mod(Gtmp + nG/2,nG)-nG/2)*ΔG

# Define the potential in direct space and construct it in momentum representation
Vext = exp.(im*x*ΔG) + exp.(-im*x*ΔG) # ≡ cos(x ΔG)
VextG = FFT(Vext)
Vin = potential_op(VextG)

# Solve the Schrodinger equation first without effective e-e interactions
eigs, vs = eigen(kinetic(0) + Vin)

bandst, wfk = erwin(Vin)

# Guess the Fermi level
EFermi = 0.0
if mod(Ne,2) == 1
  fband = trunc(Int,(Ne+1)/2)
  EFermi = bandst[trunc(Int,nk/4),fband] #nk should be a multiple of 4
else
  fband = trunc(Int,Ne/2)
  EFermi = 0.5(maximum(bandst[:,fband]) + minimum(bandst[:,fband+1]))
end
println("Fermi level: $EFermi")
println("$fband $(nk/4)")

# Calculate the population of eigenstates
occ = zeros(nk,nB)

for ik in 1:nk
  for jn in 1:nB
    occ[ik,jn] = 2occupation(bandst[ik,jn])
  end
end

@assert isapprox(sum(occ), Nes, atol=1e-3)

# Calculate the density
ρ = calc_ρ(occ,wfk)
@show ρ

using Plots
#plot(kvec,bandst[:,1:2],xlabel="k",ylabel="Energy")

#plot(x,ρ, xlabel="x/(rs bohr)",ylabel="rho")

# Coulomb potential

#plot(x,coul_pot(x,rs),xlabel="x/(rs bohr)",ylabel="v")

# Correlation energy

Δrs = 0.1
rs_array = Δrs:Δrs:40 # 39.8
vc_array = diff(Ecorr.(rs_array .- 0.05))/Δrs # creates N-1 differences
rs_array = rs_array[1:end-1] # N -> N-1
vc = CubicSplineInterpolation(rs_array,vc_array)

#p = plot(rs_array,Ecorr.(rs_array),label="onsite Corr. E")
#plot!(p,rs_array,vc.(rs_array),label="interpolated Vcorr")

# Exchange energy

rs_array = Δrs:Δrs:39.8
Δx = 0.001
xint = Δx:Δx:25
ex_vec = zeros(length(rs_array))

for (irs,rs_tmp) in enumerate(rs_array)
  ex_vec[irs] = -0.5*sum(0.5*sinc.(kF*xint/π).^2 .* coul_pot.(xint,rs_tmp)) * Δx
end

vx_array = diff(ex_vec)/Δrs # creates N-1 differences
rs_array = rs_array[1:end-1] # N -> N-1
vx = CubicSplineInterpolation(rs_array,vx_array)
ex = CubicSplineInterpolation(rs_array,ex_vec[1:end-1])

#p = plot(rs_array,ex_vec[1:end-1],label="ex")
#plot!(rs_array,vx.(rs_array),label="vx")

#plot(rs_array,Ecorr.(rs_array)./ex_vec[1:end-1],label="ecorr/ex")

# SCF

Vxc = @. ex(rs_ρ(ρ)) - vx(rs_ρ(ρ))* rs_ρ(ρ) + Ecorr(rs_ρ(ρ)) - vc(rs_ρ(ρ)) * rs_ρ(ρ)

#p = plot(Vxc,label="Vxc")
#plot!(p,real.(Vext),label="Vext")

#plot(VH(ρ) + Vxc, label="Vee") # FIXME

Vin = potential_op(VextG + FFT(VH(ρ) + Vxc))

ρ_old = ρ

residual = [1.0]
nmax = 15
mix = 0.4
ρ_all = zeros(nmax+1,nG)
ρ_all[1,:] = ρ

i = 1
while residual[end] > 5e-4 && i <= nmax
  bandst, wfk = erwin(Vin)
  global ρ = ρ*(1-mix) + mix*(calc_ρ(occ,wfk))
  Vxc = @. ex(rs_ρ(ρ)) + vx(rs_ρ(ρ))*ρ + Ecorr(rs_ρ(ρ)) + vc(rs_ρ(ρ)) * ρ
  global Vin = potential_op(VextG + FFT(VH(ρ) + Vxc))
  global i += 1
  push!(residual,sum(abs.(ρ-ρ_all[i-1,:]))/nG)
  ρ_all[i,:] = ρ
end


#p = plot(ρ,label="rho")
#plot!(p,ρ_old,label="rho_old")

#plot(ρ_all')

# Another potential

#Vext=6*[1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1]

#p = plot(Vext, label="Vext")
#plot!(p,imag.(VextG), label="Im(VextG)")
#plot!(p,3ρ, label="3rho")
