
using SciAlgs: ISBN, isvalid

@testset "ISBN validation" begin

   isbn = "84-206-8186-5"     |> ISBN
   @test isvalid(isbn)

   isbn = "84 206 8186 5"     |> ISBN
   @test isvalid(isbn)

   isbn = "978-84-473-5602-7" |> ISBN
   @test isvalid(isbn)

   isbn = "978 84 473 5602 7" |> ISBN
   @test isvalid(isbn)

   isbn = "978-84-95427-79-3" |> ISBN
   @test isvalid(isbn)

   isbn = "978 84 95427 79 3" |> ISBN
   @test isvalid(isbn)
end

