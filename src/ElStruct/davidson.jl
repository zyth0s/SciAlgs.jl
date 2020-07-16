
# Taken from https://nbviewer.jupyter.org/github/mfherbst/course_julia_day/blob/master/06_Linear_Algebra_Precision_Profiling.ipynb

using LinearAlgebra

function davidson(A, SS::AbstractArray; maxiter=100, prec=I,
                  tol=20size(A,2)*eps(eltype(A)),
                  maxsubspace=8size(SS, 2))
    m = size(SS, 2)
    for i in 1:maxiter
        Ass = A * SS

        # Use eigen specialised for Hermitian matrices
        rvals, rvecs = eigen(Hermitian(SS' * Ass))
        rvals = rvals[1:m]
        rvecs = rvecs[:, 1:m]
        Ax = Ass * rvecs

        R = Ax - SS * rvecs * Diagonal(rvals)
        if norm(R) < tol
            return rvals, SS * rvecs
        end

        println(i, "  ", size(SS, 2), "  ", norm(R))

        # Use QR to orthogonalise the subspace.
        if size(SS, 2) + m > maxsubspace
            SS = typeof(R)(qr(hcat(SS * rvecs, prec * R)).Q)
        else
            SS = typeof(R)(qr(hcat(SS, prec * R)).Q)
        end
    end
    error("not converged.")
end

# Test Davidson algorithm
# ------------------------
nev = 2
A = randn(20, 20); A = A + A' + I;

# Generate two random orthogonal guess vectors
x0 = randn(size(A, 2), nev)
x0 = Array(qr(x0).Q)

# Run the problem
E_Davidson, V_Davidson = davidson(A, x0)

# Compare with a direct diagonalization
E,V = eigen(A)
isapprox(E_Davidson, E[1:nev],   atol=1e-4) || error("Inaccurate eigenvalues")
isapprox(V_Davidson, V[:,1:nev], atol=1e-4) || error("Inaccurate eigenvectors") # sometimes fails


