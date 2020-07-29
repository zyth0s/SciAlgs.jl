
import SciAlgs.NumQuad: gauss_legendre

@testset "NumQuad: Gauss-Legendre" begin
  a = BigFloat(-1); b = BigFloat(1)
  for n in 1:6
    x, w = gauss_legendre(a,b,n)
    @test isapprox(sum(w),b-a, atol=1e-15) #"$(sum(w)) ≠ $(b-a)"
    eps = 1e-5
    if n == 1
      @test isapprox(x[1], (a+b)/2) #"$(x[1]) ≠ $((a+b)/2)"
      #
      @test isapprox(w[1], b-a) #"$(w[1]) ≠ $(b-a)"
    elseif n == 2
      @test isapprox(x[1], -0.577350,atol=eps) #"$(x[1]) ≠ $(-0.577350)"
      @test isapprox(x[2], +0.577350,atol=eps) #"$(x[1]) ≠ $(+0.577350)"
      #
      @test isapprox(w[1], 1.000000,atol=eps) #"$(w[1]) ≠ $(1.000000)"
      @test isapprox(w[2], 1.000000,atol=eps) #"$(w[1]) ≠ $(1.000000)"
    elseif n == 3
      @test isapprox(x[1], -0.774597,atol=eps) #"$(x[1]) ≠ $(-0.774597)"
      @test isapprox(x[2],  0.000000,atol=eps) #"$(x[1]) ≠ $( 0.000000)"
      @test isapprox(x[3], +0.774597,atol=eps) #"$(x[1]) ≠ $(+0.774597)"
      #
      @test isapprox(w[1], 0.555556,atol=eps) #"$(w[1]) ≠ $(0.555556)"
      @test isapprox(w[2], 0.888889,atol=eps) #"$(w[1]) ≠ $(0.888889)"
      @test isapprox(w[3], 0.555556,atol=eps) #"$(w[1]) ≠ $(0.555556)"
    elseif n == 4
      @test isapprox(x[1], -0.861136,atol=eps) #"$(x[1]) ≠ $(-0.861136)"
      @test isapprox(x[2], -0.339981,atol=eps) #"$(x[1]) ≠ $(-0.339981)"
      @test isapprox(x[3], +0.339981,atol=eps) #"$(x[1]) ≠ $(+0.339981)"
      @test isapprox(x[4], +0.861136,atol=eps) #"$(x[1]) ≠ $(+0.861136)"
      #
      @test isapprox(w[1], 0.347855,atol=eps) #"$(w[1]) ≠ $(0.347855)"
      @test isapprox(w[2], 0.652145,atol=eps) #"$(w[1]) ≠ $(0.652145)"
      @test isapprox(w[3], 0.652145,atol=eps) #"$(w[1]) ≠ $(0.652145)"
      @test isapprox(w[4], 0.347855,atol=eps) #"$(w[1]) ≠ $(0.347855)"
    elseif n == 5
      @test isapprox(x[1],-0.906180,atol=eps) #"$(x[1]) ≠ $(-0.906180)"
      @test isapprox(x[2],-0.538469,atol=eps) #"$(x[1]) ≠ $(-0.538469)"
      @test isapprox(x[3], 0.000000,atol=eps) #"$(x[1]) ≠ $( 0.000000)"
      @test isapprox(x[4],+0.538469,atol=eps) #"$(x[1]) ≠ $(+0.538469)"
      @test isapprox(x[5],+0.906180,atol=eps) #"$(x[1]) ≠ $(+0.906180)"
      #
      @test isapprox(w[1],0.236927,atol=eps) #"$(w[1]) ≠ $(0.236927)"
      @test isapprox(w[2],0.478629,atol=eps) #"$(w[1]) ≠ $(0.478629)"
      @test isapprox(w[3],0.568889,atol=eps) #"$(w[1]) ≠ $(0.568889)"
      @test isapprox(w[4],0.478629,atol=eps) #"$(w[1]) ≠ $(0.478629)"
      @test isapprox(w[5],0.236927,atol=eps) #"$(w[1]) ≠ $(0.236927)"
    elseif n == 6
      @test isapprox(x[1],-0.932470,atol=eps) #"$(x[1]) ≠ $(-0.932470)"
      @test isapprox(x[2],-0.661209,atol=eps) #"$(x[1]) ≠ $(-0.661209)"
      @test isapprox(x[3],-0.238619,atol=eps) #"$(x[1]) ≠ $(-0.238619)"
      @test isapprox(x[4],+0.238619,atol=eps) #"$(x[1]) ≠ $(+0.238619)"
      @test isapprox(x[5],+0.661209,atol=eps) #"$(x[1]) ≠ $(+0.661209)"
      @test isapprox(x[6],+0.932470,atol=eps) #"$(x[1]) ≠ $(+0.932470)"
      #
      @test isapprox(w[1],0.1713244924,atol=eps) #"$(w[1]) ≠ $(0.1713244924)"
      @test isapprox(w[2],0.360761753 ,atol=eps) #"$(w[1]) ≠ $(0.360761753 )"
      @test isapprox(w[3],0.467914    ,atol=eps) #"$(w[1]) ≠ $(0.467914    )"
      @test isapprox(w[4],0.467914    ,atol=eps) #"$(w[1]) ≠ $(0.467914    )"
      @test isapprox(w[5],0.360761753 ,atol=eps) #"$(w[1]) ≠ $(0.360761753 )"
      @test isapprox(w[6],0.1713244924,atol=eps) #"$(w[1]) ≠ $(0.1713244924)"
    end
  end
end
