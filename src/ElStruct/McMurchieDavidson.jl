
# References:
# - (HJO) "Molecular Electronic Structure Theory" 1st ed. by T. Helgaker, P. Jorgensen, J. Olsen
#         is the default reference
# - (FV) "Fundamentals of Molecular Integrals Evaluation" by J. T. Fermann, E. Valeev
# Note: This implementation is influenced by Joshua Goings Python implementation

module McMurchieDavidson

include("basis-set.jl")
import .BasisSet: BasisFunction, normalize_basis!

import LinearAlgebra: norm
import Combinatorics: binomial
import HypergeometricFunctions: drummond1F1

function E(i,j,t,Qx,a,b)
   # Hermite aexpnsion coefficients
                              # §9.5.1
                              # Note: Xpa = μ*Qx/a; Xpb = μ*Qx/b
   p = a + b                  # (9.2.10)
   μ = a*b/p                  # (9.2.12)
   if t < 0 || t > (i + j)
      return 0.0              # (9.5.5)
   elseif i == j == t == 0
      exp(-μ*Qx^2)            # (9.5.8) → (9.2.15)
   elseif j == 0
      return E(i-1,j,t-1,Qx,a,b)/(2p) -
             E(i-1,j,t,Qx,a,b)*μ*Qx/a +
             E(i-1,j,t+1,Qx,a,b)*(t+1) # (9.5.6)
   else
       return E(i,j-1,t-1,Qx,a,b)/(2p) + 
              E(i,j-1,t,Qx,a,b)*μ*Qx/b + 
              E(i,j-1,t+1,Qx,a,b)*(t+1) # (9.5.7)
   end
end

function overlap(a,lmn1,A,b,lmn2,B)
   l1,m1,n1 = lmn1
   l2,m2,n2 = lmn2
   S1 = E(l1,l2,0,A[1]-B[1],a,b)
   S2 = E(m1,m2,0,A[2]-B[2],a,b)
   S3 = E(n1,n2,0,A[3]-B[3],a,b)
   S1*S2*S3*(π/(a+b))^1.5         # (9.5.41)
end

function S(a,b)
   s = 0.0
   for (anorm,ca,aexp) in zip(a.norm,a.coefs,a.exps)
      for (bnorm,cb,bexp) in zip(b.norm,b.coefs,b.exps)
         s += anorm*bnorm*ca*cb*overlap(aexp,a.shell,a.origin,
                                        bexp,b.shell,b.origin)
      end
   end
   s
end

function kinetic(a,lmn1,A,b,lmn2,B)
   l1,m1,n1 = lmn1
   l2,m2,n2 = lmn2
   # replace (9.3.31) → (9.3.37) and group by (i,j), (i+2,j), (i-2,j)
   # i → l1, j → l2, k → m1, l → m2, m → n1, n → n2, a → b (HJO → code)
   term1 = b*(2(l2+m2+n2)+3) * overlap(a,lmn1,A,b,lmn2,B)       # (i,  j) & i,j → k,l → m,n
   term2 = -2*b^2*(        overlap(a,lmn1,A,b,(l2+2,m2,n2),B) + 
                           overlap(a,lmn1,A,b,(l2,m2+2,n2),B) +
                           overlap(a,lmn1,A,b,(l2,m2,n2+2),B) ) # (i+2,j) & i,j → k,l → m,n
   term3 = -0.5*(l2*(l2-1)*overlap(a,lmn1,A,b,(l2-2,m2,n2),B) +
                 m2*(m2-1)*overlap(a,lmn1,A,b,(l2,m2-2,n2),B) +
                 n2*(n2-1)*overlap(a,lmn1,A,b,(l2,m2,n2-2),B) ) # (i-2,j) & i,j → k,l → m,n
   term1 + term2 + term3
end

function T(a,b)
   t = 0.0
   for (ca,anorm,aexp) in zip(a.coefs,a.norm,a.exps), # a
       (cb,bnorm,bexp) in zip(b.coefs,b.norm,b.exps)  # b
      t += anorm*bnorm*ca*cb*kinetic(aexp,a.shell,a.origin, # a
                                     bexp,b.shell,b.origin) # b
   end
   t
end

function R(t,u,v,n,p,PCx,PCy,PCz,RPC)
   # Auxiliary Hermite Coulomb integrals (9.9.13)
   T = p * RPC^2
   val = 0.0
   if t == u == v == 0
      val += (-2p)^n*F_boys(n,T)                       # (9.9.14)
   elseif t == u == 0
      if v > 1
         val += (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC) # (9.9.20) 1st
      end
      val += PCz*R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)      # (9.9.20) 2nd
   elseif t == 0
      if u > 1
         val += (u-1)*R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC) # (9.9.19) 1st
      end
      val += PCy*R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)      # (9.9.19) 2nd
   else
      if t > 1
         val += (t-1)*R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC) # (9.9.18)
      end
      val += PCy*R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)      # (9.9.18)
   end
   val
end

function F_boys(n,T)
   drummond1F1(n+0.5,n+1.5,-T)/(2n+1)   # (9.8.39) → (9.7.14) 
end

function centre_of_charge(a,A,b,B)
   (a*A+b*B)/(a+b)                 # (9.2.13), (9.2.11)
end

function nuclear_attraction(a,lmn1,A,b,lmn2,B,C)
   l1,m1,n1 = lmn1
   l2,m2,n2 = lmn2
   p = a + b                     # (9.9.28)
   P = centre_of_charge(a,A,b,B) # (9.9.29)
   RPC = norm(P-C)
   
   val = 0.0
   for t in 0:l1+l2,
       u in 0:m1+m2,
       v in 0:n1+n2
      val += E(l1,l2,t,A[1]-B[1],a,b) * 
             E(m1,m2,u,A[2]-B[2],a,b) * 
             E(n1,n2,v,A[3]-B[3],a,b) * 
             R(t,u,v,0,p,P[1]-C[1],     # x
                         P[2]-C[2],     # y
                         P[3]-C[3],RPC) # z
   end
   val *= 2π/p # (9.9.40)
end

function V(a,b,C)
   v = 0.0
   for (ca,anorm,aexp) in zip(a.coefs,a.norm,a.exps), # a
       (cb,bnorm,bexp) in zip(b.coefs,b.norm,b.exps)  # b
      v += anorm*bnorm*ca*cb*nuclear_attraction(aexp,a.shell,a.origin, # a
                                                bexp,b.shell,b.origin, # b
                                                C)
   end
   v
end

function electron_repulsion(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D)
   l1,m1,n1 = lmn1
   l2,m2,n2 = lmn2
   l3,m3,n3 = lmn3
   l4,m4,n4 = lmn4
   p = a + b     # (9.9.28)
   q = c + d     # (9.9.30)
   α = p*q/(p+q) # (9.7.22)
   P = centre_of_charge(a,A,b,B) # (9.9.29)
   Q = centre_of_charge(c,C,d,D) # (9.9.31)
   RPQ = norm(P-Q)
   val = 0.0
   for t in 0:l1+l2, u in 0:m1+m2, v in 0:n1+n2,
       τ in 0:l3+l4, ν in 0:m3+m3, ϕ in 0:n3+n4
      val += E(l1,l2,t,A[1]-B[1],a,b) *
             E(m1,m2,u,A[2]-B[2],a,b) *
             E(n1,n2,v,A[3]-B[3],a,b) *
             E(l3,l4,τ,C[1]-D[1],c,d) *
             E(m3,m4,ν,C[2]-D[2],c,d) *
             E(n3,n4,ϕ,C[3]-D[3],c,d) *
             (-1)^(τ+ν+ϕ) *
             R(t+τ,u+ν,v+ϕ,0,α,P[1]-Q[1],     # x
                               P[2]-Q[2],     # y
                               P[3]-Q[3],RPQ) # z
   end
   val *= 2π^2.5/(p*q*√(p+q)) # naive (9.9.33)
end

function ERI(a,b,c,d)
   eri = 0.0
   for (ca,anorm,aexp) in zip(a.coefs,a.norm,a.exps), # a
       (cb,bnorm,bexp) in zip(b.coefs,b.norm,b.exps), # b
       (cc,cnorm,cexp) in zip(c.coefs,c.norm,c.exps), # c
       (cd,dnorm,dexp) in zip(d.coefs,d.norm,d.exps)  # d
      eri += anorm*bnorm*cnorm*dnorm*ca*cb*cc*cd*
             electron_repulsion(aexp,a.shell,a.origin, # a
                                bexp,b.shell,b.origin, # b
                                cexp,c.shell,c.origin, # c
                                dexp,d.shell,d.origin) # d
   end
   eri
end

"""
Iterator over unique indices of symmetric 4-tensor (e.g. two-electron repulsion integrals)
"""
struct twoe_iterator
   nbasis::Int
end

# https://physics.stackexchange.com/questions/87467/number-of-2-electron-integrals
Base.length(iter::twoe_iterator) = 3binomial(iter.nbasis,4) + 6binomial(iter.nbasis,3) + 4binomial(iter.nbasis,2) + binomial(iter.nbasis,1)

""" Filter non-equivalent indices by symmetry from cartesian product """
function Base.iterate(iter::twoe_iterator, state=((0, 0), (1, 1), (1, 1), (1, 1)))

   it = Iterators.product(1:iter.nbasis, 1:iter.nbasis, 1:iter.nbasis, 1:iter.nbasis)
   next = iterate(it, state)
   while next != nothing
      (i,j,k,l), state = next
      if j ≤ i && l ≤ k
         ij = i*(i+1)÷2 + j
         kl = k*(k+1)÷2 + l
         if ij >= kl
            return (i,j,k,l), state
         end
      end
      next = iterate(it, state)
   end
   nothing
end

function _twoe_iterator(bfs)
   # https://physics.stackexchange.com/questions/87467/number-of-2-electron-integrals
   nbasis = length(bfs)
   nelem = 3binomial(nbasis,4) + 6binomial(nbasis,3) + 4binomial(nbasis,2) + binomial(nbasis,1)
   list = Array{NTuple{4,Int64},1}(undef,nelem)
   ind = 1
   for i in eachindex(bfs), j in 1:i
      ij = (i*(i+1)÷2 + j)
      for k in eachindex(bfs), l in 1:k
         kl = (k*(k+1)÷2 + l)
         if ij >= kl
            list[ind] = (i,j,k,l)
            ind += 1
         end
      end
   end
   Iterators.Stateful(list)
end


end # module
