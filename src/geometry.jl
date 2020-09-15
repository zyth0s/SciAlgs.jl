
using LinearAlgebra: norm

function _cart2spherical(point::Vector{T}) where T <: Number

   length(point) == 3 || error("Must be 3D vector")
   epsrad = 1e-7
   cosθ = 0.0
   sinθ = 0.0
   cosϕ = 0.0
   sinϕ = 0.0
   #rxy = point[1]*point[1] + point[2]*point[2]
   #rxy = hypot(point[1], point[2])^2
   rxy = sum(abs2.(point[1:2]))
   if rxy < epsrad
      if point[3] >= 0.0
         cosθ = +1.0
      else
         cosθ = -1
      end
   else
      rxy = sqrt(rxy)
      cosθ = point[3] / norm(point)
      sinθ = sqrt( (1.0 - cosθ) * (1.0 + cosθ))
      cosϕ = point[1] / rxy
      sinϕ = point[2] / rxy
   end
   sinθ, cosθ, sinϕ, cosϕ
end

function cart2spherical(point::Vector{T}) where T <: Number

   length(point) == 3 || error("Must be 3D vector")

   x, y, z = point
   r = norm(point)
   ϕ = acos(z)
   #fact = hypot(x,y)
   fact = sqrt(x * x + y * y)
   if 0 < fact
      ang_x = acos(x / fact)
   else
      ang_x = acos(x)
   end
   if y < 0
      ang_x = -ang_x
   end
   θ = ang_x
   θ, ϕ, r
end


function cart2spherical(points::Matrix{T}) where T <: Number
   size(points, 2) == 3 || error("Must be Npoints × 3 matrix")
   npoints = size(points, 1)
   r = zeros(npoints)
   θ = zeros(npoints)
   ϕ = zeros(npoints)
   for i in 1:npoints
      θ[i], ϕ[i], r[i] = cart2spherical(points[i,:])
   end
   θ, ϕ, r
end


# Naive algorithm
#
points = rand(50,3)

function closest_point(points)
   npoints = size(points,1)
   minDist = Inf
   closestPair = []
   for i in 1:npoints-1
      for j in (i+1):npoints
         p = points[i,:]
         q = points[j,:]
         if norm(p-q) < minDist
            minDist = norm(p-q)
            closestPair = [p,q]
         end
      end
   end
   closestPair
end
