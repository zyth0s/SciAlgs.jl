
function read_xyz_file(filename;units="Å")
   # Read a molecule from a file with XYZ format.
   # Units are in Å.
   file = open(filename)
   natoms = parse(Int64,readline(file))
   readline(file) # Title
   atlist = []
   xyz = zeros(natoms,3)
   for i in 1:natoms
      ln = split(readline(file)) 
      push!(atlist,ln[1])
      xyz[i,:] = parse.(Float64,ln[2:4])
   end
   close(file)
   if units == "bohr"
      map(x -> x/0.52917721092,xyz)
   elseif units != "Å"
      error("Wrong specified units")
   end
   atlist, xyz
end
