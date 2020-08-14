#!/usr/bin/env julia
# Forked from https://github.com/SamChill/hartree-fock
# Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory by Szabo-Ostlund
# Appendix A and §3.5

using Formatting: printfmt
using SpecialFunctions
using LinearAlgebra

abstract type BasisFunction end

#1s Gaussian-type function
struct Gaussian1s <: BasisFunction
    #Gaussian orbital exponent
    α::Float64
    #nuclear coordinates
    center::Float64
end

#Slater Type Orbital fit with N primative gausians (STO-NG) type basis
struct STONG <: BasisFunction
    n::Integer
    #contraction coeffiecents
    d::Array{Float64} 
    #primative gaussians
    g::Array{Gaussian1s}
end

#STO-3G basis for hydrogen
function sto3g_hydrogen(center)
    sto3g(center, 1.24)
end

#STO-3G basis for helium
function sto3g_helium(center)
    sto3g(center, 2.0925)
end

#Builds a STO-3G basis that best approximates a single slater type
#orbital with Slater orbital exponent ζ
function sto3g(center, ζ)
    scaling = ζ^2
    STONG(3,[0.444635, 0.535328, 0.154329], 
            [Gaussian1s(scaling*.109818, center),
             Gaussian1s(scaling*.405771, center),
             Gaussian1s(scaling*2.22766, center)])
end

#The overlap integrals describe how the basis functions overlap
#as the atom centered gaussian basis functions are non-orthognal
#they have a non-zero overlap. The integral has the following form:
#S_{ij} = ∫ ϕ_i(r-R_a) ϕ_j(r-R_b) \mathrm{d}r
function overlap_integral(b1::STONG, b2::STONG)
    two_center_contraction(b1, b2, overlap_integral)
end

#This function calculates the overlap integral for 1s Gaussian
#orbitals using a closed-form expression.
function overlap_integral(g1::Gaussian1s, g2::Gaussian1s)
    α = g1.α
    β = g2.α
    Ra = g1.center
    Rb = g2.center

    #normalization constant
    n = (2α/π)^(3/4) * (2β/π)^(3/4)

    S  = n * (π/(α+β))^(3/2) 
    S *= exp(-α*β/(α+β) * abs(Ra-Rb)^2)

    return S
end

#
function nuclear_attraction_integral(Zc::Int, Rc::Float64, b1::STONG, 
                                     b2::STONG)
    integral(g1, g2) = nuclear_attraction_integral(Zc, Rc, g1, g2)
    two_center_contraction(b1, b2, integral)
end

function nuclear_attraction_integral(Zc::Int, Rc::Float64, g1::Gaussian1s, g2::Gaussian1s)
    α = g1.α
    β  = g2.α
    Ra = g1.center
    Rb = g2.center
    Rp = (α*Ra + β*Rb)/(α + β)

    n = (2α/π)^(3/4) * (2β/π)^(3/4)
    matrix_element  = n*-2π/(α+β)*Zc
    matrix_element *= exp(-α*β/(α+β)*abs(Ra-Rb)^2)

    t = (α+β)*abs(Rp-Rc)^2
    if abs(t) < 1e-8
        return matrix_element
    end

    matrix_element *= 0.5sqrt(π/t) * erf(sqrt(t))
    return matrix_element
end

function kinetic_energy_integral(b1::STONG, b2::STONG)
    two_center_contraction(b1, b2, kinetic_energy_integral)
end

function kinetic_energy_integral(g1::Gaussian1s, g2::Gaussian1s)
    α = g1.α
    β = g2.α
    Ra = g1.center
    Rb = g2.center

    n = (2α/π)^(3/4) * (2β/π)^(3/4)

    matrix_element  = n * α*β/(α+β)
    matrix_element *= (3-2α*β/((α+β)/abs(Ra-Rb)^2 )) 
    matrix_element *= (π/(α+β))^(3/2)
    matrix_element *= exp(-α*β/(α+β) * abs(Ra-Rb)^2)
end

function two_electron_integral(g1::STONG, g2::STONG, g3::STONG, g4::STONG)
    four_center_contraction(g1, g2, g3, g4, two_electron_integral)
end

function two_electron_integral(g1::Gaussian1s, g2::Gaussian1s, g3::Gaussian1s, 
                               g4::Gaussian1s)
    α  = g1.α
    β  = g2.α
    γ  = g3.α
    δ  = g4.α
    Ra = g1.center
    Rb = g2.center
    Rc = g3.center
    Rd = g4.center
    Rp = (α*Ra + β*Rb)/(α + β)
    Rq = (γ*Rc + δ*Rd)/(γ + δ)

    n  = (2α/π)^(3/4) * (2β/π)^(3/4)
    n *= (2γ/π)^(3/4) * (2δ/π)^(3/4)

    matrix_element  = 2n*π^(5/2)
    matrix_element /= ((α+β)*(γ+δ)*sqrt(α+β+γ+δ))
    matrix_element *= exp(-α*β/(α+β)*abs(Ra-Rb)^2 - γ*δ/(γ+δ)*abs(Rc-Rd)^2)
    t = (α+β)*(γ+δ)/(α+β+γ+δ)*abs(Rp-Rq)^2
    if abs(t) < 1e-8
        return matrix_element
    end

    matrix_element *= 0.5sqrt(π/t) * erf(sqrt(t))
    return matrix_element
end

function two_center_contraction(b1::STONG, b2::STONG, integral::Function)
    total = 0.0
    for p = 1:b1.n, q = 1:b2.n
        d1 = b1.d[p]
        d2 = b2.d[q]
        total += d1*d2*integral(b1.g[p], b2.g[q]) 
    end
    return total
end

function four_center_contraction(b1::STONG, b2::STONG, b3::STONG, b4::STONG,
                                 integral::Function)
    total = 0.0
    for p in 1:b1.n, q in 1:b2.n, r in 1:b3.n, s in 1:b4.n
        dp = b1.d[p]
        dq = b2.d[q]
        dr = b3.d[r]
        ds = b4.d[s]
        total += dp*dq*dr*ds*integral(b1.g[p], b2.g[q], b3.g[r], b4.g[s])
    end
    total
 end


function hartree_fock(R, Z)
    #println("constructing basis set")
    #ϕ = Array(BasisFunction, length(Z))
    #ϕ = Array{BasisFunction}(length(Z))
    ϕ = Array{BasisFunction}(undef,length(Z))
    for A = 1:length(Z)
        if Z[A] == 1
            ϕ[A] = sto3g_hydrogen(R[A])
        elseif Z[A] == 2
            ϕ[A] = sto3g_helium(R[A])
        end
    end
    #calculate the overlap matrix S
    #the matrix should be symmetric with diagonal entries equal to one
    #println("building overlap matrix")
    S = Matrix{Float64}(I(length(ϕ))) # eye
    for i = 1:length(ϕ)
        for j = (i+1):length(ϕ)
            S[i,j] = S[j,i] = overlap_integral(ϕ[i], ϕ[j])
        end
    end

    #calculate the kinetic energy matrix T
    #println("building kinetic energy matrix")
    T = zeros(length(ϕ), length(ϕ))
    for i = 1:length(ϕ)
        for j = i:length(ϕ)
            T[i,j] = T[j,i] = kinetic_energy_integral(ϕ[i], ϕ[j])
        end
    end

    #calculate nuclear attraction matrices V_i
    #println("building nuclear attraction matrices")
    V = zeros(length(Z), length(Z))
    for A = 1:length(Z)
        for i = 1:length(ϕ)
            for j = i:length(ϕ)
                v = nuclear_attraction_integral(Z[A], R[A], ϕ[i], ϕ[j])
                V[i,j] += v
                if i != j
                    V[j,i] += v
                end
            end
        end
    end

    #build core-Hamiltonian matrix
    #println("building core-Hamiltonian matrix")
    Hcore = T + V

    #diagonalize overlap matrix to get transformation matrix X
    #println("diagonalizing overlap matrix")
    s, U = eigen(S)
    #println("building transformation matrix")
    #X = U*diagm( 0 => s.^(-1/2))*U'
    X = U * Diagonal(s.^(-1/2)) * U'

    #calculate all of the two-electron integrals
    K = length(ϕ)
    two_electron = zeros(K,K,K,K)
    #for (μ, ν, λ, σ) in Iterators.product(1:K,1:K,1:K,1:K)
    for μ in 1:K, ν in 1:K, λ in 1:K, σ in 1:K
        coulomb  = two_electron_integral(ϕ[μ], ϕ[ν], ϕ[σ], ϕ[λ])
        two_electron[μ,ν,σ,λ] = coulomb
        exchange = two_electron_integral(ϕ[μ], ϕ[λ], ϕ[σ], ϕ[ν])
        two_electron[μ,λ,σ,ν] = exchange
    end

    P = zeros(K,K)

    total_energy = 0.0
    old_energy = 0.0
    electronic_energy = 0.0

    printfmt("{:4s} {:13s} ΔE\n", "iter", "total energy")
    for scf_iter = 1:100
        #calculate the two electron part of the Fock matrix
        G = zeros(size(Hcore))
        K = length(ϕ)

        for μ = 1:K, ν = 1:K, λ = 1:K, σ = 1:K
            coulomb  = two_electron[μ,ν,σ,λ]
            exchange = two_electron[μ,λ,σ,ν]
            G[μ,ν] += P[λ,σ]*(coulomb - 0.5exchange)
        end

        F = Hcore + G

        nuclear_energy = 0.0
        for A = 1:length(Z)
            for B = (A+1):length(Z)
                nuclear_energy += Z[A]*Z[B]/abs(R[A]-R[B])
            end
        end

        electronic_energy = 0.5dot(P,Hcore+F)
        total_energy = electronic_energy + nuclear_energy
        printfmt("{:3d} {:12.8f} {:12.4e}\n", scf_iter, total_energy, 
               total_energy - old_energy)

        if scf_iter > 2 && abs(old_energy - total_energy) < 1e-6
            break 
        end

        Fprime = X' * F * X
        epsilon, Cprime = eigen(Fprime)
        C = real(X*Cprime)

        P = 2C[:,1] * C[:,1]' # projector over occupied space

        old_energy = total_energy
    end

    return total_energy, electronic_energy
end

function heh_pes()
    file = open("../outputs/heh_pes.dat", "w")
    energy_he, energy_he = hartree_fock([0.0], [2])
    energy = []
    for r in range(0.7, 3.5, length=25)
        total_energy, electronic_energy = hartree_fock([0., r], [2, 1])
        write(file, "$r $(total_energy - energy_he)\n")
        append!(energy,total_energy)
    end
    close(file)
    plot(range(0.7,3.5,length=25),energy, xaxis="Distance [Bohr]",yaxis="Energy [Ha]")
end

if isinteractive()
   using Plots
   heh_pes()
end

