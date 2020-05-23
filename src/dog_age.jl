
# Estimate the equivalent human age of a dog


# Example
yrs = 10 # dog age
weight = 25 # kg

naive_age = 7yrs 
println("Human age estimated is $naive_age years. [naive formula]")

# https://doi.org/10.1101/829192
age1 = 16log(yrs) + 31
println("Human age estimated is $age1 years. [epigenetic]")


# https://pets.webmd.com/dogs/how-to-calculate-your-dogs-age
#print("Weight of the dog? " ); weight = parse(Int,readline(stdin))
age2 = 0

if weight < 9.071847 # 20 lbs
   if yrs == 1
      age2 = 15
   elseif yrs == 2
      age2 = 24
   elseif yrs == 3
      age2 = 28
   elseif yrs == 4
      age2 = 32
   elseif yrs == 5
      age2 = 36
   elseif yrs == 6
      age2 = 40
   elseif yrs == 7
      age2 = 44
   elseif yrs == 8
      age2 = 48
   elseif yrs == 9
      age2 = 52
   elseif yrs == 10
      age2 = 56
   elseif yrs == 11
      age2 = 60
   elseif yrs == 12
      age2 = 64
   elseif yrs == 13
      age2 = 68
   elseif yrs == 14
      age2 = 72
   elseif yrs == 15
      age2 = 76
   elseif yrs == 16
      age2 = 80
   else
      print("Error") 
   end

elseif weight < 22.67962 # 50 lbs
   if yrs == 1
      age2 = 15
   elseif yrs == 2
      age2 = 24
   elseif yrs == 3
      age2 = 28
   elseif yrs == 4
      age2 = 32
   elseif yrs == 5
      age2 = 36
   elseif yrs == 6
      age2 = 42
   elseif yrs == 7
      age2 = 47
   elseif yrs == 8
      age2 = 51
   elseif yrs == 9
      age2 = 56
   elseif yrs == 10
      age2 = 60
   elseif yrs == 11
      age2 = 65
   elseif yrs == 12
      age2 = 69
   elseif yrs == 13
      age2 = 74
   elseif yrs == 14
      age2 = 78
   elseif yrs == 15
      age2 = 83
   elseif yrs == 16
      age2 = 87
   else
      print("Error") 
   end

else
   if yrs == 1
      age2 = 15
   elseif yrs == 2
      age2 = 24
   elseif yrs == 3
      age2 = 28
   elseif yrs == 4
      age2 = 32
   elseif yrs == 5
      age2 = 36
   elseif yrs == 6
      age2 = 45
   elseif yrs == 7
      age2 = 50
   elseif yrs == 8
      age2 = 55
   elseif yrs == 9
      age2 = 61
   elseif yrs == 10
      age2 = 66
   elseif yrs == 11
      age2 = 72
   elseif yrs == 12
      age2 = 77
   elseif yrs == 13
      age2 = 82
   elseif yrs == 14
      age2 = 88
   elseif yrs == 15
      age2 = 93
   elseif yrs == 16
      age2 = 120
   else
      print("Error") 
   end
end
println("Human age estimated is $age2 years. [phenotype]")

