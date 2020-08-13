import LinearAlgebra: norm
#***********************************************************************         
#*                                                                               
#* xyzint works out the internal coordinates of a molecule.                      
#*        the "rules" for the connectivity are as follows:                       
#*        atom i is defined as being at a distance from the nearest              
#*        atom j, atom j already having been defined.                            
#*        atom i makes an angle with atom j and the atom k, which has            
#*        already been defined, and is the nearest atom to j                     
#*        atom i makes a dihedral angle with atoms j, k, and l. l having         
#*        been defined and is the nearest atom to k                              
#*                                                                               
#*        note that geo and xyz must not be the same in the call.                
#*                                                                               
#*   on input xyz    = cartesian array of numat atoms                            
#*                                                                               
#***********************************************************************         
function xyzint(xyz)

   numat = size(xyz,2)
   na = zeros(Int,size(xyz)[2])
   nb = zeros(Int,size(xyz)[2])
   nc = zeros(Int,size(xyz)[2])
   geo = zeros(3,size(xyz)[2])

   #nai1 = 0
   #nai2 = 0
   k = 0
   for i in 1:numat
      na[i] = 2
      nb[i] = 3
      nc[i] = 4
      im1 = i - 1
      if im1 == 0 
         @goto noatoms
      end
      sum = 100.0
      for j in 1:im1
         R = (xyz[1,i] - xyz[1,j])^2 +
             (xyz[2,i] - xyz[2,j])^2 +
             (xyz[3,i] - xyz[3,j])^2
         R = norm(xyz[:,i]-xyz[:,j])
         if R < sum && na[j] != j && nb[j] != j
            sum = R
            k = j
         end
      end
      # i is nearest to k
      na[i] = k
      if i > 2
         nb[i] = na[k]
      end
      if i > 3 
         nc[i] = nb[k]
      end
      @label noatoms
   end
   na[1] = 0
   nb[1] = 0
   nc[1] = 0
   nb[2] = 0
   nc[2] = 0
   nc[3] = 0
   geo = xyzgeo(xyz,numat,na,nb,nc)
   for i in 1:numat
      if geo[3,i] > 180.0
         geo[3,i] -= 360.0
      end
   end
   geo, na, nb, nc
end

#***********************************************************************         
#*                                                                               
#*   xyzgeo converts coordinates from cartesian to internal.                     
#*                                                                               
#*     on input xyz  = array of cartesian coordinates                            
#*              numat= number of atoms                                           
#*              na   = numbers of atom to which atoms are related                
#*                     by distance                                               
#*              nb   = numbers of atom to which atoms are related                
#*                     by angle                                                  
#*              nc   = numbers of atom to which atoms are related                
#*                     by dihedral                                               
#*                                                                               
#*    on output geo  = internal coordinates in angstroms, radians,               
#*                     and radians                                               
#*                                                                               
#***********************************************************************         
function xyzgeo(xyz,numat,na,nb,nc)
   geo = zeros(3, numat)
   for i in 2:numat
      j = copy(na[i])
      k = copy(nb[i])
      l = copy(nc[i])
      if i < 3
         @goto l1
      end
      ii = i
      geo[2,i] = bangle(xyz,ii,j,k)
      geo[2,i] = rad2deg(geo[2,i])
      if i < 4
         @goto l1
      end
      geo[3,i] = dihedral(xyz,ii,j,k,l)
      geo[3,i] = rad2deg(geo[3,i])
      @label l1
      geo[1,i] = norm(xyz[:,i] - xyz[:,j])
   end
   geo[1,1] = 0.0
   geo[2,1] = 0.0
   geo[3,1] = 0.0
   geo[2,2] = 0.0
   geo[3,2] = 0.0
   geo[3,3] = 0.0
   geo
end

#*********************************************************************           
#*                                                                               
#* bangle calculates the angle between atoms i,j, and k. The
#*        cartesian coordinates are in xyz.                                      
#*                                                                               
#*********************************************************************           
function bangle(xyz,i,j,k)

   d2ij = (xyz[1,i]-xyz[1,j])^2+
          (xyz[2,i]-xyz[2,j])^2+
          (xyz[3,i]-xyz[3,j])^2
   d2jk = (xyz[1,j]-xyz[1,k])^2+
          (xyz[2,j]-xyz[2,k])^2+
          (xyz[3,j]-xyz[3,k])^2
   d2ik = (xyz[1,i]-xyz[1,k])^2+
          (xyz[2,i]-xyz[2,k])^2+
          (xyz[3,i]-xyz[3,k])^2
   xy = sqrt(d2ij*d2jk)
   temp = 0.5(d2ij+d2jk-d2ik) / xy
   if temp >  1.0
      temp=1.0
   end
   if temp < -1.0
      temp=-1.0
   end
   angle = acos( temp )
end

#*********************************************************************
#*                                                                               
#*      dihedral calculates the dihedral angle between atoms i, j, k,  
#*            and l.  the cartesian coordinates of these atoms        
#*            are in array xyz.                                      
#*                                                                  
#*     dihedral is a modified version of a subroutine of the same name
#*           which was written by dr. w. theil in 1973.              
#*                                                                    
#*********************************************************************
function dihedral(xyz,i,j,k,l)

      xi1 = xyz[1,i] - xyz[1,k]                                                     
      xj1 = xyz[1,j] - xyz[1,k]                                                     
      xl1 = xyz[1,l] - xyz[1,k]                                                     
      yi1 = xyz[2,i] - xyz[2,k]                                                     
      yj1 = xyz[2,j] - xyz[2,k]                                                     
      yl1 = xyz[2,l] - xyz[2,k]                                                     
      zi1 = xyz[3,i] - xyz[3,k]                                                     
      zj1 = xyz[3,j] - xyz[3,k]                                                     
      zl1 = xyz[3,l] - xyz[3,k]                                                     
      # rotate around z axis to put kj along y axis                              
      dist= sqrt(xj1^2 + yj1^2 + zj1^2)                                          
      cosa=zj1/dist                                                             
      if cosa > 1.0
         cosa=1.0                                              
      end
      if cosa < -1.0
         cosa=-1.0                                            
      end
      ddd=1.0-cosa^2                                                         
      if ddd <= 0.0
         @goto l10                                                   
      end
      yxdist=dist* sqrt(ddd)                                                    
      if yxdist > 1.0e-9
         @goto l20                                             
      end
   @label l10
      xi2=xi1                                                                   
      xl2=xl1                                                                   
      yi2=yi1                                                                   
      yl2=yl1                                                                   
      costh=cosa                                                                
      sinth=0.0                                                                
      @goto l30                                                                  
      @label l20
      cosph=yj1/yxdist                                                          
      sinph=xj1/yxdist                                                          
      xi2=xi1*cosph-yi1*sinph                                                   
      xj2=xj1*cosph-yj1*sinph                                                   
      xl2=xl1*cosph-yl1*sinph                                                   
      yi2=xi1*sinph+yi1*cosph                                                   
      yj2=xj1*sinph+yj1*cosph                                                   
      yl2=xl1*sinph+yl1*cosph                                                   
#      rotate kj around the x axis so kj lies along the z axis                  
      costh=cosa                                                                
      sinth=yj2/dist                                                            
      @label l30
      yi3=yi2*costh-zi1*sinth                                                   
      yl3=yl2*costh-zl1*sinth                                                   
      angle = dangle(xl2,yl3,xi2,yi3)
      if angle < 0.0
         angle+=2π
      end
      if angle >= 2π 
         angle=0.0                                   
      end
      angle
end

#**********************************************************************          
#*                                                                               
#*    dangle  determines the angle between the points (a1,a2), (0,0),             
#*          and (b1,b2).  the result is put in rcos.                             
#*                                                                               
#**********************************************************************          
function dangle(a1, a2, b1, b2)

   if abs(a1) < 1e-6 && abs(a2) < 1e-6
      return 0.0                        
   end
   if abs(b1) < 1e-6 && abs(b2) < 1e-6
      return 0.0                        
   end
   anorm=1.0/ sqrt(a1^2 + a2^2)                                            
   bnorm=1.0/ sqrt(b1^2 + b2^2)                                            
   a1 = a1 * anorm
   a2 = a2 * anorm
   b1 = b1 * bnorm
   b2 = b2 * bnorm
   sinth = (a1*b2) - (a2*b1)                                                     
   costh = a1*b1 + a2*b2                                                         
   if costh > 1.0
      costh = 1.0
   end
   if costh < -1.0
      costh = -1.0
   end
   rcos= acos(costh)                                                         
   if  abs(rcos) < 4.0e-4
      return 0.0                                         
   end
   if sinth > 0.0
      rcos = 6.2831853 - rcos
   end
   rcos = -rcos
end

# Example
if isinteractive()
   xyz = zeros(3,8)
   xyz[:,1] = [ 0.77389  -0.01480 -0.00796]
   xyz[:,2] = [-0.77389   0.01480 -0.00796]
   xyz[:,3] = [-1.18197  -0.85433 -0.51417]
   xyz[:,4] = [-1.14844   0.89890 -0.51417]
   xyz[:,5] = [-1.16521   0.02228  1.00446]
   xyz[:,6] = [ 1.18197   0.85433 -0.51417]
   xyz[:,7] = [ 1.14844  -0.89890 -0.51417]
   xyz[:,8] = [ 1.16521  -0.02228  1.00446]
   elnames = ["C","C","H","H","H","H","H","H"]

   geo, na, nb, nc = xyzint(xyz)

   zmat = hcat(elnames,na,geo[1,:],nb,geo[2,:],nc,geo[3,:])
end
