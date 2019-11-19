#!/usr/bin/env julia

using LinearAlgebra
using Formatting: printfmt


# Metric tensor (Å²)
function metric_tensor(a,b,c,α,β,γ)
   G = [ a*a        a*b*cos(γ) a*c*cos(β);
         a*b*cos(γ) b*b        b*c*cos(α);
         a*c*cos(β) b*c*cos(α) c*c       ]
end

function characterize_metric_tensor(G)
   V = volume_unit_cell(G) # (Å³)
   println(" * Metric tensor G:")
   printfmt("   {:7.3f} {:7.3f} {:7.3f} \n",G[1,1],G[1,2],G[1,3])
   printfmt("   {:7.3f} {:7.3f} {:7.3f} \n",G[2,1],G[2,2],G[2,3])
   printfmt("   {:7.3f} {:7.3f} {:7.3f} \n",G[3,1],G[3,2],G[3,3])
   printfmt(" * Volume = {:15.4f} aₒ³ \n",V)
end

function unit_volume_unit_cell(a,b,c,α,β,γ)
   # α,β,γ in radians
   sqrt(1 + 2*cos(α)*cos(β)*cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2)
end
function volume_unit_cell(a,b,c,α,β,γ)
   # α,β,γ in radians
   a*b*c*sqrt(1 + 2*cos(α)*cos(β)*cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2)
end

function volume_unit_cell(G)
   sqrt(det(G)) # (Å³)
end

function lattice_vectors(a,b,c,α,β,γ)
   # α,β,γ in radians
   unitV = unit_volume_unit_cell(a,b,c,α,β,γ)
   aᵣ = sin(α) / (a * unitV)
   cos_γᵣ = (cos(α)*cos(β) - cos(γ))/ (sin(α)*sin(β))
   sin_γᵣ = sqrt(1.0 - cos_γᵣ^2)
   a1 = [ 1.0/ aᵣ; -cos_γᵣ / sin_γᵣ / aᵣ; cos(β) * a]
   a2 = [       0;            b * sin(α); b * cos(α)]
   a3 = [       0;                     0;          c]
   a1, a2, a3
end

function recip_cell_parameters(a,b,c,α,β,γ)
   V = volume_unit_cell(a,b,c,α,β,γ)
   aᵣ = b * c * sin(α)/V
   bᵣ = a * c * sin(β)/V
   cᵣ = a * b * sin(γ)/V
   αᵣ = acos((cos(β)*cos(γ)-cos(α))/(sin(β)*sin(γ)))
   βᵣ = acos((cos(α)*cos(γ)-cos(β))/(sin(γ)*sin(α)))
   γᵣ = acos((cos(α)*cos(β)-cos(γ))/(sin(β)*sin(α)))
   aᵣ,bᵣ,cᵣ,αᵣ,βᵣ,γᵣ
end

function reciprocal_vectors(a1,a2,a3)
   Tr = 2*π*inv([a1 a2 a3])'
   b1 = Tr[:,1]
   b2 = Tr[:,2]
   b3 = Tr[:,3]
   b1, b2, b3
end

function reciprocal_vectors(T)
   Tr = 2*π*inv(T)'
end

normcryst(r,G) = sqrt(r'*G*r)
innercryst(r1,r2,G) = r1'*G*r2
anglecryst(r1,r2,G) = acos(innercryst(r1,r2,G)/(normcryst(r1,G)*normcryst(r2,G)))

function interplanar_spacing(h,k,l, Gr)
   # d[hkl] = 1 / || rᵣ || = 1 / sqrt[ (h k l) G* (h k l)T ]
   d = [h; k; l]' * Gr * [h; k; l]
   sqrt(1.0 / d)
end

function interplanar_spacing(h,k,l,aᵣ,bᵣ,cᵣ,αᵣ,βᵣ,γᵣ)
   abc2 = h^2*aᵣ^2 + k^2*bᵣ^2 + l^2*cᵣ^2
   bc = h*k*cos(γᵣ)
   ac = h*l*cos(βᵣ)                                  
   ba = k*l*cos(αᵣ)                                  
   d2 = abc2 + 2.0*(aᵣ*bᵣ*bc + aᵣ*cᵣ*ac + bᵣ*cᵣ*ba)          
   sqrt(1.0/d2)                                    
end



function cart_frac(a,b,c,α,β,γ)
   V = volume_unit_cell(a,b,c,α,β,γ)
   # Cartesian -> fractional coordinates x' = M x
   M = [ 1/a  -cos(γ)/(a*sin(γ))  b*c*(cos(α)*cos(γ)-cos(β))/(V*sin(γ));
           0        1/(b*sin(γ))  a*c*(cos(β)*cos(γ)-cos(α))/(V*sin(γ));
           0                   0                           a*b*sin(γ)/V]
   # Fractional -> cartesian coordinates x = Minv x'
   Minv = [ a b*cos(γ) c*cos(β)                       ;
            0 b*sin(γ) c*(cos(α)-cos(β)*cos(γ))/sin(γ);
            0 0        V/(a*b*sin(γ))                 ]
   M, Minv
end

function map2Seitz(mapping)
   println(" * Mapping: ", mapping...)
   w = [float(subs(c, (x,0.), (y,0.), (z,0.))) for c in mapping]
   W = zeros(3,3)
   W[:,1] = [float(subs(c, (x,1.), (y,0.), (z,0.))) for c in mapping] - w
   W[:,2] = [float(subs(c, (x,0.), (y,1.), (z,0.))) for c in mapping] - w
   W[:,3] = [float(subs(c, (x,0.), (y,0.), (z,1.))) for c in mapping] - w
   W, w
end

function characterize_Seitz(W,w,G)
   println(" * Seitz symbol (W|w):")
   printfmt("   {:7.3f} {:7.3f} {:7.3f} | {:7.3f} \n",W[1,1],W[1,2],W[1,3],w[1])
   printfmt("   {:7.3f} {:7.3f} {:7.3f} | {:7.3f} \n",W[2,1],W[2,2],W[2,3],w[2])
   printfmt("   {:7.3f} {:7.3f} {:7.3f} | {:7.3f} \n",W[3,1],W[3,2],W[3,3],w[3])
   W'*G*W ≈ G || error("Not an isometry")
   detW = det(W)
   traceW = tr(W)
   if detW ≈ 1.0
      printfmt("   det(W) = {:+7.3f}: Rotation \n",detW)
      if traceW != 3
         v, u = eigen(W)
         for i in 1:3
            if v[i] == 1
               u = u[:,i]
            end
         end
      end
      if traceW == 3
         type = 1
      end
      if traceW == 2
         type = 6
      end
      if traceW == 1
         type = 4
      end
      if traceW == 0
         type = 3
      end
      if traceW == -1
         type = 2
      end
      k = type
   elseif detW ≈ -1
      printfmt("   det(W) = {:+7.3f}: Rotoinversion \n",detW)
      if traceW != -3
         v, u = eigen(W)
         for i in 1:3
            if v[i] == -1
               u = u[:,i]
            end
         end
      end
      if traceW == -3
         type = -1
         k = 2
      end
      if traceW == -2
         type = -6
         k = 6
      end
      if np.trace(W) == -1
         type = -4
         k = 4
      end
      if traceW == 0
         type = -3
         k = 6
      end
      if traceW == 1
         type = -2
         k = 2
      end
   else
      error("  det(W) = {:+7.3f}",detW)
   end

   W^k ≈ I || error("Order is wrong!")

   tmp = I
   for i in 1:k
      tmp += W^i
   end
   screw_glide_component = tmp * w / k

   printfmt("   Crystallographic symmetry operation: {:2d} (order {:2d})\n",type,k)
   if type == 1
      if w ≈ zeros(3)
         printfmt("    Identity\n")
      else
         printfmt("   Translation vector {:7.3f} {:7.3f} {:7.3f}\n",w...)
      end
   end
   if type == 6
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
      if !isapprox(screw_glide_component,zeros(3,3))
         printfmt("   Screw component:  {:7.3f} {:7.3f} {:7.3f}\n",screw_glide_component...)
      end
   end
   if type == 4
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
      if !isapprox(screw_glide_component,zeros(3,3))
         printfmt("   Screw component:  {:7.3f} {:7.3f} {:7.3f}\n",screw_glide_component...)
      end
   end
   if type == 3
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
      if !isapprox(screw_glide_component,zeros(3,3))
         printfmt("   Screw component:  {:7.3f} {:7.3f} {:7.3f}\n",screw_glide_component...)
      end
   end
   if type == 2
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
   end
   if type == -1
      printfmt("   Point of inversion:  {:7.3f} {:7.3f} {:7.3f}\n".format(0.5*w[0],0.5*w[1],0.5*w[2]))
   end
   if type == -6
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
   end
   if type == -4
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
   end
   if type == -3
      printfmt("   Fixed axis: {:7.3f} {:7.3f} {:7.3f}\n",u...)
   end
   if type == -2
      printfmt("   Reflection plane perpendicular to {:7.3f} {:7.3f} {:7.3f}\n",u...)
      if !isapprox(screw_glide_component,zeros(3,3))
         printfmt("   Glide component:  {:7.3f} {:7.3f} {:7.3f}\n",screw_glide_component...)
      end
   end
end

function transform_point(r,W,w,M,Minv, normalize=false)
   # p is a point in cartesian coordinates
   p = M*p
   if normalize
      p[:] -= floor(p[:])
   end
   p = W*p + w
   Minv*p
end

# ================================================================================

#function example1()
# Example 1: 
a = 13.90222379    # Å
b = 8.08681840     # Å     
c = 16.0706089     # Å
α = 90 * π / 180  # radians
β = 90 * π / 180  # radians
γ = 90 * π / 180  # radians

G = metric_tensor(a,b,c,α,β,γ)
characterize_metric_tensor(G)


M, Minv = cart_frac(a,b,c,α,β,γ)

# Cartesian to crystallographic
x = [3.23629868,2.0217046, -6.16388206] # column vector of atom 
M*x

using SymPy
@vars x y z

mapping = -x + 0.5, -y, z + 0.5
W, w = map2Seitz(mapping)
characterize_Seitz(W,w,G)

a1, a2, a3 = lattice_vectors(a,b,c,α,β,γ)
T = [a1 a2 a3]
println(" * Lattice vectors (a1 a2 a3) = T:")
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",a1[1],a2[1],a3[1])
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",a1[2],a2[2],a3[2])
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",a1[3],a2[3],a3[3])
b1, b2, b3 = reciprocal_vectors(a1,a2,a3)
Tr = [b1 b2 b3]
println(" * Reciprocal vectors (b1 b2 b3) = T*:")
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",Tr[1,1],Tr[1,2],Tr[1,3])
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",Tr[2,1],Tr[2,2],Tr[2,3])
printfmt("   {:7.3f} {:7.3f} {:7.3f}\n",Tr[3,1],Tr[3,2],Tr[3,3])

#end

#function example2()
#   # Example 2: FeO₂...
#   a = 9.691    # Å
#   b = 8.993     # Å     
#   c = 5.231     # Å
#   α = 90     * π / 180  # radians
#   β = 108.61 * π / 180  # radians
#   γ = 90     * π / 180  # radians
#
#   r1 = [ 0.7431, 0.5154, 0.2770] # fractional
#   r2 = [ 0.6313, 0.5162, -0.1178] # fractional
#   r3 = [ 0.8679, 0.3378,  0.1812] # fractional
#
#   G = metric_tensor(a,b,c,α,β,γ)
#
#   # Distance between two fractional coordinates
#   dist = normcryst(r2-r1,G) # Å
#   println("Distance Fe-O is ",dist)
#   dist = normcryst(r3-r1,G) # Å
#   println("Distance Fe-O2 is ",dist)
#
#   # Angle between two vectors in cryst coordinates
#   angulo = anglecryst(r2-r1,r3-r1,G)
#   println("Angulo O-Fe-O2 is ",angulo * 180 / π)
#end

#example1()
#example2()

# TODO(dmc)
# full     Hermann-Maugin
# short    Hermann-Maugin
# extended Hermann-Maugin

function interplanar_spacing(h,k,l,a,b,c,α,β,γ)
   unitV = unit_volume_unit_cell(a,b,c,α,β,γ)
   unitV / (h^2/a^2*sin(α)^2 
          + k^2/b^2*sin(β)^2 
          + l^2/c^2*sin(γ)^2 
          + 2*k*l/(b*c)*(cos(β)*cos(γ) - cos(α)) 
          + 2*h*l/(a*c)*(cos(γ)*cos(α)-cos(β)) 
          + 2*h*k/(a*b)*(cos(α)*cos(β)-cos(γ))    )
end

function interplanar_spacing(θ,λ)
   " Uses Bragg's Law"
   λ / (2.0 * sin(θ) )
end

function calc_θ(d,λ)
   # Bragg's Law: 2d = n λ sin(θ)
   θ = 0.0
   if 2d < λ
      return d, θ
   end
	tea = 0.5λ/d
	tea = tea/sqrt(1.0-tea*tea)
	θ = atan(tea)
end


function crystallite_size(β,λ,θ,K=0.9)
   " L = K * λ / (β * (2*θ) * cos(θ) )"
   " K is the Scherrer constant"
   " β is the peak width in degrees (FWHM)"
   " θ is the peak position in degrees"
   " λ is the wavelength of X-ray radiation in Å"
   " Most effective for L < 200 nm"
   " Return nanometers"
   K * λ / (β * (2*θ) * cosd(θ) )
end
