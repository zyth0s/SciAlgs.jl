
# Calculates interreticular distances excluding
# triclinic systems.
#
# Inspired by ERADIS, a code made accessible by  
# A. Le Bail - March 1998. lebail@univ-lemans.fr
#
# Daniel Menendez Crespo @ MPI CPFS, Dresden

using Formatting: printfmt


#Codes according to the Bravais lattice and crystal system
bravais_codes = Dict(
               "CubicP"          =>  [1,2,3,0,0,0,0],
               "CubicI"          =>  [1,2,3,0,0,0,7],
               "CubicF"          =>  [1,2,3,4,5,6,0],
               "TetragonalP"     =>  [1,0,3,0,0,0,0],
               "TetragonalI"     =>  [1,0,3,0,0,0,7],
               "HexagonalP"      =>  [1,0,3,0,0,0,0],
               "OrthorhombicP"   =>  [1,0,0,0,0,0,0],
               "OrthorhombicI"   =>  [1,0,0,0,0,0,7],
               "OrthorhombicF"   =>  [1,0,0,4,5,6,0],
               "OrthorhombicA"   =>  [1,0,0,4,0,0,0],
               "OrthorhombicB"   =>  [1,0,0,0,5,0,0],
               "OrthorhombicC"   =>  [1,0,0,0,0,6,0],
               "RhomboedricR"    =>  [0,2,3,0,0,0,0],
               "MonoclinicP"     =>  [0,0,0,0,0,0,0],
               "MonoclinicA"     =>  [0,0,0,4,0,0,0],
               "MonoclinicC"     =>  [0,0,0,0,0,6,0])

function volume_unit_cell(a,b,c,α,β,γ)
   # α,β,γ in radians
   a*b*c*sqrt(1 + 2*cos(α)*cos(β)*cos(γ) - cos(α)^2 - cos(β)^2 - cos(γ)^2)
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

function metric_tensor(a,b,c,α,β,γ)
   G = [ a*a        a*b*cos(γ) a*c*cos(β);
         a*b*cos(γ) b*b        b*c*cos(α);
         a*c*cos(β) b*c*cos(α) c*c       ]
end

function interplanar_spacing(h,k,l, Gr)
   # d[hkl] = 1 / || rᵣ || = 1 / sqrt[ (h k l) G* (h k l)T ]
   d = [h; k; l]' * Gr * [h; k; l]
   sqrt(1.0 / d)
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

function find_reflections(spg_info, params)
   a, b, c, α, β, γ, bravais_lattice = spg_info
   dmin, dmax, hmax, kmax, lmax, λ = params

   aᵣ,bᵣ,cᵣ,αᵣ,βᵣ,γᵣ = recip_cell_parameters(a,b,c,α,β,γ)
   Gr = metric_tensor(aᵣ,bᵣ,cᵣ,αᵣ,βᵣ,γᵣ)
   bravais_lattice == "Triclinic" && 
                  @warn("The program does not work for triclinic cells")
   bravais_code = bravais_codes[bravais_lattice]
   h = 0
   k = 0
   l = 0
   nel = 0
   num = 0
   hkl_table = zeros(Int64,600,3)
   dis = zeros(600)
   θs = zeros(600)
   # Welcome to GOTO HELL!! Spaguetti code ahead ...
   @label label10; 
   i = 1
   @label label20; ip = bravais_code[i]
   if ip == 0
      @goto label8
   elseif ip == 1
      @goto label1
   elseif ip == 2
      @goto label2
   elseif ip == 3
      @goto label3
   elseif ip == 4
      @goto label4
   elseif ip == 5
      @goto label5
   elseif ip == 5
      @goto label5
   elseif ip == 6
      @goto label6
   elseif ip == 6
      @goto label6
   elseif ip == 7
      @goto label7
   end

   @label label1; if l < 0 
      @goto label12
   end
   @goto label8
   @label label2; if k < abs(l)
      @goto label12
   end
   @goto label8
   @label label3; if h < k
      @goto label12
   end
   @goto label8
   @label label4; mh = k+l
   @goto label171
   @label label5; mh = h+l
   @goto label171
   @label label6; mh = h + k
   @goto label171
   @label label7; mh = h + k + l
   @label label171; if mod(mh,2) != 0 
      @goto label12
   end
   @label label8; i += 1
   if i > 7
      @goto label9
   end
   @goto label20
   @label label9 
   d = interplanar_spacing(h,k,l,Gr)
   θ = calc_θ(d,λ)
   if d < dmin 
      @goto label12
   end
   if d > dmax
      @goto label12
   end
   num += 1
   hkl_table[num,:] = [h,k,l]
   dis[num] = d
   θs[num] = θ
   if num == 600
      @goto label18
   end
   @label label12; if l == 0
      @goto label13
   end
   if h == 0
      @goto label13
   end
   l = -l
   if nel == 1
      @goto label13
   end
   nel = 1
   @goto label10
   @label label13; 
   nel = 0
   l += 1
   if l <= lmax
      @goto label10
   end
   l = 0
   k += 1
   if k <= kmax
      @goto label10
   end
   l = 0
   k = 0
   h += 1
   if h <= hmax
      @goto label10
   end
   @label label18
   return num, dis, θs, hkl_table
end

function show_reflections(num,dis,θs,hkl_table)
   println(" Number of calculated planes = $num")
   n = num - 1
   if (n - 1) < 0
      return
   end

   sortinds = sortperm(dis, rev=true)
   hkl_table = hkl_table[sortinds,:]
   θs = θs[sortinds]
   dis = dis[sortinds]

   println("---- ---- ---- ---------- -------")
   println("   h    k    l       d[Å]   2θ[°]")
   println("---- ---- ---- ---------- -------")
   for j in 1:num
      h, k, l, d, θ = hkl_table[j,:]..., dis[j], rad2deg(θs[j])
      printfmt("{:4d} {:4d} {:4d} {:10.5f} {:7.3f}\n",h, k, l, d, 2θ)
   end
   println("---- ---- ---- ---------- -------")
end
 
#-----------------------------------------------------
# Example 1: CaCuF4 tetragonal I
a = 5.3770    # Å
b = 5.3770    # Å     
c = 10.3200   # Å
α = deg2rad(90) # radians
β = deg2rad(90) # radians
γ = deg2rad(90) # radians
bravais_lattice = "TetragonalI"

dmin = 1.20 # Å
dmax = 12.0 # Å
hmax = 19
kmax = 19
lmax = 19
λ = 1.54056 # Cu Kα (Å)

#-----------------------------------------------------

spg_info = [a, b, c, α, β, γ, bravais_lattice]
params = [dmin, dmax, hmax, kmax, lmax, λ]
num, dis, θs, hkl_table = find_reflections(spg_info, params)
show_reflections(num, dis, θs, hkl_table)

