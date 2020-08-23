
import HTTP
import JSON

function read_xyz_file(filename;units="Å")
   # Read a molecule from a file with XYZ format.
   # Units are in Å.
   file = open(filename)
   natoms = parse(Int64,readline(file)) # 1st line
   readline(file)                       # 2nd line: Title
   atlist = String[]
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

function fetch_basis(basis_name,atlist)
   r = HTTP.request("GET","http://basissetexchange.org/api/basis/" * basis_name * "/format/json",query=Dict("elements" => atlist))
   r.status == 200 || error("Could not obtain data from the BSE. Check the error information above")
   JSON.parse(String(r.body))
end
