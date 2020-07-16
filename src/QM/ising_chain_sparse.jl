#=
Quantum Ising model on chain of N spins 

   H = -∑_⟨i,j⟩ σᵢᶻ ⊗ σⱼᶻ - h ∑ᵢ σᵢˣ 

* subindex indicates site
* superindex x,z indicates which Pauli matrix
* h is the amplitude of the transverse magnetic field in x direction
* ⟨i,j⟩ only nearest neighbors interactions
* σᶻ|↓⟩ = |↑⟩
* σˣ|↑⟩ = |↓⟩

Tutorial given by Carsten Bauer, Katharine Hyatt
https://juliaphysics.github.io/PhysicsTutorials.jl/tutorials/general/quantum_ising/quantum_ising.html
=#

using LinearAlgebra
using SparseArrays
using ArnoldiMethod

const ⊗ = kron

function eigen_sparse(x)
    decomp, history = partialschur(x, nev=1, which=SR()); # only solve for the ground state
    vals, vecs = partialeigen(decomp);
    return vals, vecs
end

#=
Binary `BitArray` representation of the given integer `num`, padded to length `N`.
=#
bit_rep(num::Integer, N::Integer) = BitArray(parse(Bool, i) for i in string(num, base=2, pad=N))

#=
    generate_basis(N::Integer) -> basis

Generates a basis (`Vector{BitArray}`) spanning the Hilbert space of `N` spins.
=#
function generate_basis(N::Integer)
    nstates = 2^N
    basis = Vector{BitArray{1}}(undef, nstates)
    for i in 0:nstates-1
        basis[i+1] = bit_rep(i, N)
    end
    return basis
end

function TransverseFieldIsing_sparse(;N,h)
    @assert N < 21 "The chain is too long"
    id = [1 0; 0 1] |> sparse
    σˣ = [0 1; 1 0] |> sparse
    σᶻ = [1 0; 0 -1] |> sparse

    first_term_ops = fill(id, N)
    first_term_ops[1] = σᶻ
    first_term_ops[2] = σᶻ

    second_term_ops = fill(id, N)
    second_term_ops[1] = σˣ

    H = spzeros(Int, 2^N, 2^N) # note the spzeros instead of zeros here
    for i in 1:N-1
        H -= foldl(⊗, first_term_ops)
        first_term_ops = circshift(first_term_ops,1)
    end

    for i in 1:N
        H -= h*foldl(⊗, second_term_ops)
        second_term_ops = circshift(second_term_ops,1)
    end
    H
end

function magnetization(state)
    N = Int(log2(length(state)))
    M = 0.
    for i in 1:length(state)
        bstate = bit_rep(i-1,N)
        bstate_M = 0.
        for spin in bstate
            bstate_M += (state[i]^2 * (spin ? 1 : -1))/N
        end
        @assert abs(bstate_M) <= 1
        M += abs(bstate_M)
    end
    return M
end

# Test

basis = generate_basis(8)
@time H = TransverseFieldIsing_sparse(N=8, h=100)
vals, vecs = eigen_sparse(H)

groundstate = vecs[:,1];
abs2.(groundstate)

magnetization(groundstate,basis) ≈ 0.2748106008973601

using PyPlot

hs = 10 .^ range(-2., stop=2., length=10)
Ns = 2:2:20
xlabel("h")
ylabel("M(h)")
@time for N in Ns
    M = zeros(length(hs))
    for (i,h) in enumerate(hs)
        H = TransverseFieldIsing_sparse(N=N, h=h)
        vals, vecs = eigen_sparse(H)
        groundstate = @view vecs[:,1]
        M[i] = magnetization(groundstate)
    end
    semilogx(hs,M,label="N = $N")
    println(M)
end
legend(loc="lower left")

