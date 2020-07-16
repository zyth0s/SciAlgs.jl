#=
Quantum Ising model on chain of N spins 

   H = -∑_{i,j}^3 σᵢˣ ⊗ σⱼˣ - h ∑ᵢ^3 σᵢᶻ 

* subindex indicates site
* superindex x,z indicates which Pauli matrix
* h is the amplitude of the transverse magnetic field in z direction
* ⟨i,j⟩ only nearest neighbors interactions
* σᶻ|↓⟩ = |↑⟩
* σˣ|↑⟩ = |↓⟩

Tutorial given by Guifre Vidal
https://github.com/eschnett/2018-computational-physics-course/blob/master/exact-diagonalization-module/2018-computational-physics-course-ExactDiag1.ipynb
=#

using LinearAlgebra

const ⊗ = kron

id = [1 0; 0 1]
σˣ = [0 1; 1  0] # Pauli matrices
σᶻ = [1 0; 0 -1]

#       1    2    3      1    2    3      1    2    3
HXX =  σˣ ⊗ σˣ ⊗ id  +  id ⊗ σˣ ⊗ σˣ  +  σˣ ⊗ id ⊗ σˣ  
HZ  =  σᶻ ⊗ id ⊗ id  +  id ⊗ σᶻ ⊗ id  +  id ⊗ id ⊗ σᶻ

H_max = 2 # maximum magnetic field
nsteps = 20
h = range(0,stop=H_max,length=nsteps)

Eg = zeros(nsteps)
for i in 1:nsteps
   H = HXX + h[i]*HZ
   D,U = eigen(H)
   Eg[i] = first(D) #[1] # take the lowest
end

using PyPlot

plot(h,Eg,marker=".",color="b")
grid("on")
title("Ground state energy: quantum Ising model for N=3 spins")
xlabel("magnetic field h")
ylabel("Energy")

