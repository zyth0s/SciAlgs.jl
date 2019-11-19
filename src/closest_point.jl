
# Naive algorithm
#
using LinearAlgebra
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
