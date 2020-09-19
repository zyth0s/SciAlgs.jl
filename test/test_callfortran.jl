

@testset "Mixed language: Call Fortran" begin

   include("../src/interface_fortran/iso_c/callfortran.jl")
   @test v[] == 4

end

