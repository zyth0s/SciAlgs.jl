
using SciAlgs.NumQuad: clenshaw_curtis

@testset "NumQuad: Clenshaw-Curtis" begin
  for i in [1,2,3,4,5,7]
    for a in BigFloat[-4,-3,-2,-1, 0], b in BigFloat[1,2,3,4,5]
      x, w = clenshaw_curtis(a,b,i)
      tol = 1e-15
      @test isapprox(sum(w),b-a,atol=tol) #"$(sum(w)) ≠ $(b-a)"
      if i == 1
        @test isapprox(w[1], b-a,atol=tol) #"$(w[1]) ≠ $(b-a)"
      elseif i == 2
        @test isapprox(w[1], (b-a)/2,atol=tol) #"$(w[1]) ≠ $((b-a)/2)"
        @test isapprox(w[2], (b-a)/2,atol=tol) #"$(w[1]) ≠ $((b-a)/2)"
      elseif i == 3
        @test isapprox(w[1],   (b-a)/6,atol=tol) #"$(w[1]) ≠ $(  (b-a)/6)"
        @test isapprox(w[2], 4*(b-a)/6,atol=tol) #"$(w[1]) ≠ $(4*(b-a)/6)"
        @test isapprox(w[3],   (b-a)/6,atol=tol) #"$(w[1]) ≠ $(  (b-a)/6)"
      elseif i == 4
        @test isapprox(w[1],   (b-a)/18,atol=tol) #"$(w[1]) ≠ $(   (b-a)/18)"
        @test isapprox(w[2], 8*(b-a)/18,atol=tol) #"$(w[1]) ≠ $( 8*(b-a)/18)"
        @test isapprox(w[3], 8*(b-a)/18,atol=tol) #"$(w[1]) ≠ $( 8*(b-a)/18)"
        @test isapprox(w[4],   (b-a)/18,atol=tol) #"$(w[1]) ≠ $(   (b-a)/18)"
      elseif i == 5                                                    
        @test isapprox(w[1],   (b-a)/30,atol=tol) #"$(w[1]) ≠ $(   (b-a)/30)"
        @test isapprox(w[2], 8*(b-a)/30,atol=tol) #"$(w[1]) ≠ $( 8*(b-a)/30)"
        @test isapprox(w[3],12*(b-a)/30,atol=tol) #"$(w[1]) ≠ $(12*(b-a)/30)"
        @test isapprox(w[4], 8*(b-a)/30,atol=tol) #"$(w[1]) ≠ $( 8*(b-a)/30)"
        @test isapprox(w[5],   (b-a)/30,atol=tol) #"$(w[1]) ≠ $(   (b-a)/30)"
      elseif i == 7
        @test isapprox(w[1],  9*(b-a)/630,atol=tol) #"$(w[1]) ≠ $(  9*(b-a)/630)"
        @test isapprox(w[2], 80*(b-a)/630,atol=tol) #"$(w[1]) ≠ $( 80*(b-a)/630)"
        @test isapprox(w[3],144*(b-a)/630,atol=tol) #"$(w[1]) ≠ $(144*(b-a)/630)"
        @test isapprox(w[4],164*(b-a)/630,atol=tol) #"$(w[1]) ≠ $(164*(b-a)/630)"
        @test isapprox(w[5],144*(b-a)/630,atol=tol) #"$(w[1]) ≠ $(144*(b-a)/630)"
        @test isapprox(w[6], 80*(b-a)/630,atol=tol) #"$(w[1]) ≠ $( 80*(b-a)/630)"
        @test isapprox(w[7],  9*(b-a)/630,atol=tol) #"$(w[1]) ≠ $(  9*(b-a)/630)"
      end
    end
  end
end

