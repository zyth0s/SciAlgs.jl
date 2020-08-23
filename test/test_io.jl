
using SciAlgs: fetch_basis

@testset "I/O: Fetch basis from BSSE" begin

   basis_set = fetch_basis("sto-3g",["H"])
   @test basis_set["name"] == "STO-3G"
   @test basis_set["elements"]["1"]["electron_shells"][1]["exponents"][1] == "0.3425250914E+01"
   @test basis_set["elements"]["1"]["electron_shells"][1]["exponents"][2] == "0.6239137298E+00"
   @test basis_set["elements"]["1"]["electron_shells"][1]["exponents"][3] == "0.1688554040E+00"
end

