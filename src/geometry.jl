
using LinearAlgebra

function cart2polar(point)
  epsrad = 1e-7
  cost = 0.0
  sint = 0.0
  cosp = 0.0
  sinp = 0.0
  rxy = point[1]*point[1] + point[2]*point[2]
  if rxy < epsrad
    if point[3] >= 0.0
      cost = +1.0
    else
      cost = -1
    end
  else
    rxy = sqrt(rxy)
    cost = point[3] / norm(point)
    sint = sqrt( (1.0 - cost) * (1.0 + cost))
    cosp = point[1] / rxy
    sinp = point[2] / rxy
  end
  cost, sint, cosp, sinp
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
