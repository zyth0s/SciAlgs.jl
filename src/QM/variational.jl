
# "A linear variational exercise with a simple non-orthogonal basis for the
#  particle-in-the-box problem" Eur. J. Phys. 31 (2010) 101–114, 
# doi:10.1088/0143-0807/31/1/010
# Linear variational calculation of the particle in the box using 
# ϕₙ(x) = Nₙ * (x * (a-x))ⁿ as basis functions
#
using SpecialFunctions
using LinearAlgebra
using Printf
#using Plots
import PyPlot
pyplt = PyPlot
mpl = PyPlot.matplotlib
pyplt.matplotlib.style.reload_library()
pyplt.matplotlib.style.use("5in_color")
mpl.use(backend="Qt5Agg")

function overlap(r,s)
  # overlap integral : ⟨ϕᵣ|ϕₛ⟩
  beta(r+s+1,r+s+1) / sqrt(beta(2*r+1,2*r+1)*beta(2*s+1,2*s+1))
end

function kinetic(r,s)
  # Kinetic energy integral: ⟨ϕᵣ|T|ϕₛ⟩
  # ϵ = h² / (8ma²) is used as energy scale of the problem
  rs = r+s;
  (r*s/π^2) * (beta(rs-1,rs-1) - 4*beta(rs,rs-1) + 
     4*beta(rs+1,rs-1)) / sqrt(beta(2*r+1,2*r+1)*beta(2*s+1,2*s+1))
end

function phi(n,y)
  # Evaulate "Nₙ y (1-y)" for a grid of values of y in [0,1]
  # "a" is used as the length scale of the problem.
  # Notice the use of "element-by-element" operations
  norma = 1/sqrt(beta(2*n+1,2*n+1))
  (y .* (1 .-y)) .^ n .* norma
end

# Dimension of the basis set: 1, 2, ... MAXN
# Build up the H and S matrices:
MAXN = 3
H = zeros(MAXN,MAXN); S = zeros(MAXN,MAXN);
for i = 1:MAXN
  for j = 1:MAXN
    H[i,j] = kinetic(i,j); S[i,j] = overlap(i,j);
  end
end

# Solve the generalized eigenproblem. Use symmetric orthogonalization
V = S^(-1/2);
Hprime = V' * H * V;
W, Cprime = eigen(Hprime);
C = V * Cprime;

# Sort (ascending) the eigenvalues and reorder the C matrix accordingly:
ival = sortperm(W);
Csort = C[:,ival];
val = W[ival]
#uout = open("results.dat","w")
open("results.dat","w") do uout
  write(uout,"Sorted eigenvalues:\n")
  @printf(uout,"%12.9f %12.9f %12.9f\n",val[1:MAXN]...)
end
#close(uout)

# Evaulate the basis functions on a grid
y = range(0,stop=1,length=101);
Phi = zeros(MAXN,101);
for i = 1:MAXN
  Phi[i,:] = phi(i,y);
end

# Evaluate and plot eigenvectors
EigVect = C' * Phi;
#grid("on"); xlabel("y = x/a"); ylabel("Phi / sqrt(a)")
#plot(y, EigVect[1,:], grid=true, xlabel = "y = x/a", ylabel = "Phi / sqrt(a)")
#plot!(y, EigVect[2,:])
#plot!(y, EigVect[3,:])
#savefig("results.png")
#print("results.ps","-color");

fig = pyplt.figure()
fig.set_dpi(260)
ax = fig.add_subplot(111)
# Plot data
ax.plot(y,EigVect[1,:],label="Eigvec 1")
ax.plot(y,EigVect[2,:],label="Eigvec 2")
ax.plot(y,EigVect[3,:],label="Eigvec 3")
# Set axis labels
ax.set_xlabel("y = x/a")
ax.set_ylabel("Phi / sqrt(a)")
# Add legend to plot
ax.legend(bbox_to_anchor=(0.7, 0.3), loc="upper right", ncol=1,frameon=false)
# Save plot
fig.tight_layout(pad=0.1)
pyplt.savefig("../figures/nonortho_variational.pdf")

